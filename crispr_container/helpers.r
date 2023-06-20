bitmap_pred <- function(fasta_df, fasta_index = 1, batch_size = 200, return_summary = FALSE, step = 300,
                        char_sequence = NULL, model = NULL, rc = FALSE, range = NULL, 
                        include_label = FALSE, summary_method = "mean") {
  
  if (return_summary) stopifnot(summary_method %in% c("mean", "max"))
  
  maxlen <- model$input_shape[[2]]
  
  if (is.null(char_sequence)) {
    k <- fasta_index
    char_sequence <- fasta_df$Sequence[k]
  }
  
  char_len <- nchar(char_sequence)
  
  if (!is.null(range)) {
    stopifnot(length(range) == 2)
    stopifnot(range[1] < range[2])
    stopifnot(range[2] <= nchar(char_sequence))
    char_sequence <- substr(char_sequence, range[1], range[2])
    shift_position <- range[1] - 1
  } else {
    shift_position <- 0
  }
  
  if (nchar(char_sequence) < maxlen) {
    stop("Sequence too short")
  }
  
  start_ind <- seq(1, nchar(char_sequence) - maxlen + 1, by = step)
  char_sequence_sub <- substr(char_sequence, 
                              min(start_ind),
                              max(start_ind) + maxlen - 1)
  
  x <- seq_encoding_label(sequence = NULL,
                          maxlen = maxlen, 
                          vocabulary = c("a", "c", "g", "t"),
                          start_ind = start_ind,
                          ambiguous_nuc = "zero",
                          nuc_dist = NULL,
                          quality_vector = NULL,
                          use_coverage = FALSE,
                          max_cov = NULL,
                          cov_vector = NULL, 
                          n_gram = NULL, 
                          n_gram_stride = 1,
                          char_sequence = char_sequence_sub)
  
  if (rc) {
    x <- x[ , , 4:1]
    x <- x[ , dim(x)[2]:1, ]
  }
  
  pred_list <- list()
  num_samples <- dim(x)[1]
  num_eval_steps <- ceiling(num_samples/batch_size) 
  for (i in 1:num_eval_steps) {
    index_start <- ((i-1) * batch_size) + 1
    index_end <- min((i * batch_size), num_samples)
    index <- index_start : index_end
    pred_input <- x[index, , ]
    if (length(index) == 1) pred_input <- array(pred_input, dim = c(1, dim(pred_input)))
    
    pred_tensor <- predict(model, pred_input)
    pred_list[[i]] <- keras::array_reshape(pred_tensor, dim = c(1, prod(dim(pred_tensor))))
  }
  
  pred <- do.call(cbind, pred_list) %>% as.vector()
  
  pos_list <- list()
  for (i in 1:length(start_ind)) {
    pos_list[[i]] <- seq(start_ind[i], start_ind[i] + maxlen - 1) 
  }
  position <- unlist(pos_list) + shift_position
  
  if (rc) {
    position <- char_len - position + 1
  }
  
  df <- data.frame(pred = pred, position = position)
  
  if (include_label) {
    df$label <- rep(1:num_samples, each = maxlen)
  }
  
  if (return_summary & summary_method == "mean") {  
    df <- df %>% dplyr::group_by(position) %>% dplyr::summarise(pred = mean(pred))
  }
  
  if (return_summary & summary_method == "max") {  
    df <- df %>% dplyr::group_by(position) %>% dplyr::summarise(pred = max(pred))
  }
  
  df <- data.frame(position = df$position, pred = df$pred)
  df
  
}

# crispr_gap: After what numbber of negative predictions to start new CRISPR region.
# conf_cutoff: Confidence decision threshold.
# pos_rate: What percentage of predictions must be above conf_cutoff.
# min_seq_len: Minimum sequence length.
filter_crispr <- function(states_df,
                          crispr_gap = 10,
                          conf_cutoff = 0.5,
                          pos_rate = 0.8,
                          min_seq_len = 120,
                          maxlen = 200) {
  
  stopifnot(all(c("conf_CRISPR", "conf_non_CRISPR", "position") %in% colnames(states_df)))
  if (any(duplicated(states_df$position))) {
    stop("positions should all be unique (position column in states_df)")
  }
  step <- states_df$position[2] - states_df$position[1]
  crispr_list <- list()
  states_df <- states_df %>% dplyr::filter(conf_CRISPR > conf_cutoff)
  states_df <- states_df[order(states_df$position), ]
  row_num <- nrow(states_df)
  if (row_num == 0) {
    print("All confidence scores below conf_cutoff")
    return(NULL)
  } 
  
  crispr_index <- 1
  crispr_start <- states_df$position[1]
  
  if (row_num > 1) {
    for (i in 1:(row_num-1)) {
      
      l_name <- paste0("CRISPR_region_", crispr_index)
      current_pos <- states_df$position[i]
      next_pos <- states_df$position[i+1]
      
      if ((abs(current_pos - next_pos) > crispr_gap) & (i != (row_num-1))) {
        
        index <- (states_df$position >= crispr_start) & (states_df$position <= current_pos)
        crispr_list[[l_name]] <- states_df[index, ]
        crispr_start <- next_pos
        crispr_index <- crispr_index + 1
      }
      
      if (i == (row_num-1)) {
        if (abs(current_pos - next_pos) <= crispr_gap) {
          index <- states_df$position >= crispr_start
          crispr_list[[l_name]] <- states_df[index, ]
        } else {
          index <- (states_df$position >= crispr_start) & (states_df$position <= current_pos)
          crispr_list[[l_name]] <- states_df[index, ]
          
          # single sample at end 
          crispr_index <- crispr_index + 1
          l_name <- paste0("CRISPR_region_", crispr_index)
          crispr_list[[l_name]] <- states_df[nrow(states_df), ]
        }
      }
    } 
  } else {
    l_name <- paste0("CRISPR_region_", crispr_index)
    crispr_list[[l_name]] <- states_df
  }
  
  # filter by positivity rate
  for (i in names(crispr_list)) {
    df <- crispr_list[[i]]
    seq_len <- df$position[nrow(df)] - df$position[1]
    ## consider step size
    num_possible_pos_pred <- ((seq_len - 1)/step) + 1
    cov_rate <- nrow(df)/num_possible_pos_pred 
    if (cov_rate < pos_rate) {
      crispr_list[[i]] <- NULL
    } 
  }
  
  # filter by size 
  for (i in names(crispr_list)) {
    df <- crispr_list[[i]]
    seq_len <- df$position[nrow(df)] - df$position[1] + 1
    if (seq_len < min_seq_len) {
      crispr_list[[i]] <- NULL
    } 
  }
  
  for (i in names(crispr_list)) {
    df <- crispr_list[[i]]
    crispr_list[[i]] <- df
  }
  crispr_list
}