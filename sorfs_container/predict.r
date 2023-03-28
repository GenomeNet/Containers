#!/usr/bin/env Rscript
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3) # suppress TF messages
suppressPackageStartupMessages(suppressWarnings(library(deepG)))
library(magrittr)
library(optparse)
library(ggplot2)


option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "/data/example.fasta"),
  make_option(c("-o", "--output"), type = "character", default = "/data/output.tsv"),
  make_option(c("-s", "--step"), type = "integer", default = 1),
  make_option(c("-b", "--batch_size"), type = "integer", default = 32,
help = "Number of samples processed in one batch [default %default]",
metavar = "number"));

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

args <- commandArgs(trailingOnly = TRUE)

# Check if output directory exists and create
if (!file.exists("/data/output")) {
  dir.create("/data/output")
}

model <- keras::load_model_hdf5("/data/model_binary.h5", compile = FALSE)

bitmap_pred <- function(fasta_df = NULL, fasta_path = NULL, fasta_index = 1, batch_size = 200, 
                        return_summary = FALSE,
                        step = 1, out_names = NULL,
                        char_sequence = NULL, model = NULL, rc = FALSE, range = NULL, 
                        include_label = FALSE, summary_method = "mean") {
  
  if (!is.null(fasta_path)) fasta_df <- microseq::readFasta(fasta_path)
  if (return_summary) stopifnot(summary_method %in% c("mean", "max"))
  
  maxlen <- model$input_shape[[2]]
  num_out_layers <- length(model$outputs)
  
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
  
  x <- deepG::seq_encoding_label(sequence = NULL,
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
    l <- list()
    for (k in 1:num_out_layers) {
      l[[k]] <- keras::array_reshape(pred_tensor[[k]], 
                                     dim = c(1, prod(dim(pred_tensor[[k]]))))
      
    }
    pred_list[[i]] <- l
  }
  
  per_output_pred <- list() 
  for (j in 1:num_out_layers) {
    l <- list()
    for (i in 1:length(pred_list)) {
      l[[i]] <- pred_list[[i]][[j]]
    }
    per_output_pred[[j]] <- do.call(cbind, l) %>% as.vector()
  }
  
  pos_list <- list()
  for (i in 1:length(start_ind)) {
    pos_list[[i]] <- seq(start_ind[i], start_ind[i] + maxlen - 1) 
  }
  position <- unlist(pos_list) + shift_position
  
  if (rc) {
    position <- char_len - position + 1
  }
  
  df_pred <- do.call(cbind, per_output_pred) %>% as.data.frame()
  if (is.null(out_names)) {
    names(df_pred) <- paste0("class_", 1:num_out_layers)
  } else {
    names(df_pred) <- out_names
  }
  df <- data.frame(df_pred, position = position)
  
  if (include_label) {
    df$label <- rep(1:num_samples, each = maxlen)
  }
  
  if (return_summary) {  
    l <- list()
    for (i in 1:ncol(df_pred)) {
      col_name <- names(df_pred)[i]
      if (return_summary & summary_method == "mean") {
        df_temp <- df %>% dplyr::group_by(position) %>% dplyr::summarise(pred = mean(get(col_name)))
      } else {
        df_temp <- df %>% dplyr::group_by(position) %>% dplyr::summarise(pred = max(get(col_name)))
      }
      names(df_temp) <- c("position", col_name)
      if (i > 1) df_temp[["position"]] <- NULL
      l[[i]] <- df_temp
    }
    df <- do.call(cbind, l)
  }
  
  df
}

out_names <- c("positive_direction", "negative_direction")

df <- bitmap_pred(fasta_df = NULL,
                  fasta_path = opt$input,
                  fasta_index = 1,
                  batch_size = opt$batch_size, 
                  return_summary = FALSE,
                  step = opt$step,
                  out_names = out_names,
                  char_sequence = NULL,
                  model = model,
                  rc = FALSE,
                  range = NULL, # c(1, 2000), 
                  include_label = FALSE, 
                  summary_method = "mean")

head(df)

df2 <- data.frame(position = rep(df$position, 2),
                  conf = c(df$positive_direction, df$negative_direction),
                  direction = rep(c("positive", "negative"), each = nrow(df)))

write.table(df2, file = "/data/output/predictions.tsv", sep = "\t", quote = F, row.names = FALSE)
