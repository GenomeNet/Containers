#!/usr/bin/env Rscript
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3) # suppress TF messages
suppressPackageStartupMessages(suppressWarnings(library(deepG)))
library(magrittr)
library(optparse)
library(ggplot2)
#
suppressPackageStartupMessages(suppressWarnings(library(zoo)))
suppressPackageStartupMessages(suppressWarnings(library(plotly)))
suppressPackageStartupMessages(suppressWarnings(library(htmlwidgets)))

# Function to create a ggplot for a specific entry
create_plot <- function(data, entry_name, threshold = 0.1, min_length = 5) {
    subset <- data[which(data$name == entry_name), ]
    # Add rolling mean
    window_size <- 10  # Adjust this value for the desired window size
    subset$rolling_mean <- rollmean(subset$is_virus,
      window_size, fill = NA, align = "center")

    # Create the ggplot
    p <- ggplot() +
        # Add the line plot
        geom_line(data = subset, aes(x = element_num,
          y = rolling_mean),
          color = "black", linewidth = .5) +
        # Add the text layer
        geom_text(data = subset[1, ], aes(x = element_num, y = 1,
          label = paste("Sequence:", entry_name)),
          hjust = 0, vjust = 1) +
        scale_y_continuous(limits = c(0, 1))
    # Convert ggplot to interactive plotly plot
    ggplotly(p, tooltip = c("x", "y"))
}


option_list <- list(
  make_option(c("-m", "--model"), type = "character", default = NULL,
              help = "model to use (genus/binary)", metavar = "character"),
  make_option(c("-i", "--input"), type = "character", default = "/data/example.fasta"),
  make_option(c("-s", "--step"), type = "integer", default = 1000,
help = "Step size to iterate though sequences [default %default]",
metavar = "number"),
  make_option(c("-b", "--batch_size"), type = "integer", default = 32,
help = "Number of samples processed in one batch [default %default]",
metavar = "number"),
  make_option(c("-f", "--fast"), action = "store_true", default = FALSE,
help = "FAST mode, only predict one sample per FASTA entry"),
  make_option(c("-e", "--by_entry"), action = "store_true", default = FALSE,
              help="Toggle to make predictions by FASTA entry instead of whole file"));

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if output directory exists and create
if (!file.exists("/data/output")) {
  dir.create("/data/output")
}

a <- Sys.time()

# Load file from stdin and save to disk
message("Loading file")

args <- commandArgs(trailingOnly = TRUE)

if (!opt$model %in% c("binary", "genus")){
	message("please select the model via --model. Supported models are binary and genus")
}

if (opt$model == "binary") {
	model <- keras::load_model_hdf5("virus_binary_2023-01-23.hdf5",
    compile = FALSE)
	message("Processing file")

  if (opt$by_entry) {
  
    if (opt$fast) {
      
      h5_path <- "/data/output/predictions.h5"
      predict_model(output_format = "one_pred_per_entry",
                    model = model,
                    layer_name = "dense_3",
                    path_input = opt$input,
                    filename = h5_path,
                    batch_size = opt$batch_size,
                    output_type = "h5",
                    padding = "standard",
                    mode = "label",
                    round_digits = 4,
                    verbose = FALSE,
                    include_seq = FALSE)
      pred_list <- load_prediction(h5_path = h5_path)
     

						pred_list <- load_prediction(h5_path = h5_path)
						df <- as.data.frame(pred_list$states)
						colnames(df) <- c("virus", "non_virus")
						message(paste0("Numer of FASTA entries predicted as viral: ",length(which(df$virus >= 0.5))))
						message(paste0("Numer of FASTA entries predicted as non-viral: ",length(which(df$virus < 0.5))))               
						write.table(df, file = "/data/output/predictions.tsv", sep = "\t", quote =F, row.names = FALSE)
				#		p <- ggplot(data = df, aes(x = virus)) + geom_histogram(bins = 30)
					#	p <- p + xlab("Virus probability (%)") + ylab("Number of FASTA entries")
					#	ggplotly(p)
						# Save the interactive plot as an HTML file
					#	htmlwidgets::saveWidget(p, "/data/output/predictions_binary_histogram.html")

  }  else {

    h5_path <- "/data/output/predictions.h5"
    predict_model(output_format = "by_entry_one_file",
                  model = model,
                  layer_name = "dense_3",
                  path_input = opt$input,
                  filename = h5_path,
                  batch_size = opt$batch_size,
                  output_type = "h5",
                  padding = "standard",
                  mode = "label",
                  round_digits = 4,
                  verbose = FALSE,
                  include_seq = FALSE)
    pred_list <- load_prediction(h5_path = h5_path)
    l <- list()
    for (i in 1:length(pred_list)) {
      l[[i]] <- data.frame(entry_id = i,
                          name = names(pred_list)[[i]],
                          virus_mean = round(mean(pred_list[[i]]$states[,1]), digits = 3),
                          virus_sd = round(sd(pred_list[[i]]$states[,1]), digits = 3),
                          virus_min = round(min(pred_list[[i]]$states[,1]), digits = 3),
                          virus_max = round(max(pred_list[[i]]$states[,1]), digits = 3),
                          num_pred = length(pred_list[[i]]$states[,1]))
    }
    pred_all <- do.call("rbind", l)
    write.table(pred_all, file = "/data/output/predictions_binary.tsv", sep = "\t", quote = F, row.names = FALSE)
   message("Numer of FASTA entries predicted as viral: ", length(which(pred_all$virus_mean >= .5)))
		 message("Numer of FASTA entries predicted as non-viral: ", length(which(pred_all$virus_mean < .5)))

    l <- list()
  for (i in 1:length(pred_list)) {
    l[[i]] <- data.frame(name = names(pred_list)[[i]],
                        is_virus = pred_list[[i]]$states[,1],
                        element_num = seq_along(pred_list[[i]]$states[,1]))
  }
 # write.table(pred_all, file = "/data/output/predictions_binary_summary.tsv", sep = "\t", quote = F, row.names = FALSE)


  # Combine the dataframes in the list into a single dataframe
  pred_all <- do.call("rbind", l)
write.table(pred_all, file = "/data/output/predictions_binary_summary.tsv", sep = "\t", quote = F, row.names = FALSE)

    # make plot


    # Create a list to store individual interactive plots
    plot_list <- lapply(unique(pred_all$name), function(entry_name) {
      create_plot(pred_all, entry_name)
    })

    # Combine the plots vertically using subplot
    combined_plots <- do.call(subplot, c(plot_list, list(nrows = length(plot_list), margin = 0.05)))

    # Save the interactive plot as an HTML file
    htmlwidgets::saveWidget(combined_plots, "/data/output/predictions_binary.html")

 }
  } else {
    pred <- predict_model(output_format = "one_seq",
                        model = model,
                        layer_name = "dense_3",
                        sequence = NULL,
                        path_input = opt$input,
                        round_digits = 4,
                        filename = NULL,
                        step = opt$step,
                        batch_size = opt$batch_size,
                        verbose = FALSE,
                        return_states = TRUE,
                        padding = "standard",
                        mode = "label",
                        format = "fasta")
   # b = Sys.time()
 
    df <- data.frame(pred$states)
    names(df) <- c("non_viral", "viral")
    write.table(df, file = "/data/output/predictions_binary.tsv", sep = "\t", row.names = FALSE)
      agg <- colMeans(df)
    agg_o <- agg[order(agg, decreasing = T)]
    write.table(agg_o, file = "/data/output/predictions_binary_summary.tsv", sep = "\t", row.names = TRUE)
    agg_o <- agg_o[which(round(agg_o, digits = 2) > 0)]
    for (i in 1:length(agg_o)){
      message(paste0("Predicted FASTA file as ", names(agg_o[i]), " (" , round(agg_o[i]* 100, digits = 1), "%)"))
    }
   # message(paste0("Prediction took ", round(as.numeric(difftime(time1 = b, time2 = a, units = "secs")), 2), " seconds"))
  }
}

if (opt$model == "genus") {

	model <- keras::load_model_hdf5("virus_genus_2023-01-23.hdf5", compile = FALSE)
	genus_labels <- readRDS("virus_genus_2023-01-23_labels.rds")

 if (opt$by_entry) {
 #  message("Processing by entry")
    if (opt$fast) {
 # message("FAST mode")

		 h5_path <- "/data/output/predictions.h5"
      predict_model(output_format = "one_pred_per_entry",
                    model = model,
                    layer_name = NULL,
                    path_input = opt$input,
                    filename = h5_path,
                    batch_size = opt$batch_size,
                    output_type = "h5",
                    padding = "standard",
                    mode = "label",
                    round_digits = 4,
                    verbose = FALSE,
                    include_seq = FALSE)
      pred_list <- load_prediction(h5_path = h5_path)
    names(pred_list)
    df <- as.data.frame(pred_list$states)
    genus_labels <- readRDS("virus_genus_2023-01-23_labels.rds")
    colnames(df) <- genus_labels
    write.table(df, file = "/data/output/predictions.tsv", sep = "\t", quote = F, row.names = FALSE)


# function to get the name of the prediction with the highest value
get_max_prediction <- function(row) {
  max_value <- max(row)
  max_index <- which(row == max_value)
  message(paste0("Predicted FASTA entry as ", colnames(df)[max_index]))
  return(list(prediction = colnames(df)[max_index], value = round(max_value, 2)))
}

# apply the function to each row of the matrix
max_predictions_list <- apply(df, 1, get_max_prediction)
# Convert list to data.frame
max_predictions <- do.call(rbind, lapply(max_predictions_list, as.data.frame))
# Save max_predictions as a tsv file
write.table(max_predictions, file = "/data/output/predictions_genus_summary.csv", sep = "\t", quote = F, row.names = F)




				} else {

					 h5_path <- "/data/output/predictions.h5"
    predict_model(output_format = "by_entry_one_file",
                  model = model,
                  layer_name = "dense_3",
                  path_input = opt$input,
                  filename = h5_path,
                  batch_size = opt$batch_size,
                  output_type = "h5",
                  padding = "standard",
                  mode = "label",
                  round_digits = 4,
                  verbose = FALSE,
                  include_seq = FALSE)
   pred_list <- load_prediction(h5_path = h5_path)

      l <- list()
    for (i in 1:length(pred_list)) {
        df <- as.data.frame(pred_list[[i]]$states)
        names(df) <- genus_labels
        agg <- colMeans(df)
        agg_o <- agg[order(agg, decreasing = T)]
        message(paste0("Predicted FASTA entry as ", names(agg_o[i]), " (", round(agg_o[i], digits = 1),")" ))
        l[[i]] <- data.frame(name = names(pred_list)[[i]],
                            predicted = names(agg_o[i]),
                            probability = round(agg_o[i], digits = 3))
    }

    pred_all <- do.call("rbind", l)
    rownames(pred_all) <- NULL
  
    write.table(df, file = "/data/output/predictions.tsv", sep = "\t", row.names = FALSE)


				}
	} else {

	message("Processing file")
	pred <- predict_model(output_format = "one_seq",
                      model = model,
                      layer_name = "dense_3", 
                      sequence = NULL,
                      path_input = opt$input,
                      round_digits = 4,
                      filename = NULL,
                      step = opt$step,
                      batch_size = opt$batch_size,
                      verbose = FALSE,
                      return_states = TRUE,
                      padding = "standard",
                      mode = "label",
                      format = "fasta")

	#message("Writing predictions")
	df <- data.frame(pred$states)
	names(df) <- genus_labels
	write.table(df, file = "/data/output/predictions_genus.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
	b = Sys.time()
	
	agg <- colMeans(df)
	agg_o <- agg[order(agg, decreasing = T)]
  write.table(agg_o, file = "/data/output/predictions_genus_summary.tsv", sep = "\t", row.names = TRUE, quote = FALSE)
	message("Top 5 predictions of the sample:")
  for (i in 1:5){
    message(paste0("Predicted FASTA sample as ", names(agg_o[i]), " (" , round(agg_o[i]* 100, digits = 1), "%)"))
	}
  #message(paste0("Prediction took ", round(as.numeric(difftime(time1 = b, time2 = a, units = "secs")), 3), " seconds"))

}
}