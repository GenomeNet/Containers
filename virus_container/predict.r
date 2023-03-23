#!/usr/bin/env Rscript
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3) # suppress TF messages
suppressPackageStartupMessages(suppressWarnings(library(deepG)))
library(magrittr)
library(optparse)
library(ggplot2)

option_list = list(
  make_option(c("-m", "--model"), type="character", default = NULL,
              help = "model to use (genus/binary)", metavar = "character"),
  make_option(c("-i", "--input"), type = "character", default = "/data/example.fasta"),
  make_option(c("-s", "--step"), type = "integer", default = 1000,
help = "Step size to iterate though sequences [default %default]",
metavar = "number"),
  make_option(c("-b", "--batch_size"), type = "integer", default = 32,
help = "Number of samples processed in one batch [default %default]",
metavar = "number"),
  make_option(c("-e", "--by_entry"), action="store_true", default = FALSE,
              help="Toggle to make predictions by FASTA entry instead of whole file"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Check if output directory exists and create
if(!file.exists("/data/output")) {
  dir.create("/data/output")
}

a = Sys.time()

# Load file from stdin and save to disk
message("Loading file")

args = commandArgs(trailingOnly = TRUE)

if (!opt$model %in% c("binary", "genus")){
	message("please select the model via --model. Supported models are binary and genus")
}

if (opt$model == "binary"){
	model <- keras::load_model_hdf5("virus_binary_2023-01-23.hdf5", compile = FALSE)
	message("Processing file")

  if (opt$by_entry){
    message("Processing by entry")
    h5_path <- "/data/output/states.h5"
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
    write.table(pred_all, file = "/data/output/virus_per_entry.tsv", sep = "\t", row.names = FALSE)
  
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
    b = Sys.time()
    message("Writing predictions")
    df <- data.frame(pred$states)
    names(df) <- c("non_viral", "viral")
    write.csv(df, file = "/data/output/prediction_binary.csv", row.names = FALSE)
      agg <- colMeans(df)
    agg_o <- agg[order(agg, decreasing = T)]
    write.csv(agg_o, file = "/data/output/prediction_summary_binary.csv", row.names = TRUE)
    agg_o <- agg_o[which(round(agg_o, digits = 2) > 0)]
    for (i in 1:length(agg_o)){
      message(paste0("The sample appears to be ", names(agg_o[i]), " (" , round(agg_o[i]* 100, digits = 1), "%)"))
    }
    message(paste0("Prediction took ", round(as.numeric(difftime(time1 = b, time2 = a, units = "secs")), 2), " seconds"))
  }
}

if (opt$model == "genus"){
	model <- keras::load_model_hdf5("virus_genus_2023-01-23.hdf5", compile = FALSE)
	genus_labels <- readRDS("virus_genus_2023-01-23_labels.rds")
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

	message("Writing predictions")
	df <- data.frame(pred$states)
	names(df) <- genus_labels
	write.csv(df, file = "/data/output/prediction_genus.csv", row.names = FALSE)
	b = Sys.time()
	
	agg <- colMeans(df)
	agg_o <- agg[order(agg, decreasing = T)]
  write.csv(agg_o, file = "/data/output/prediction_summary_genus.csv", row.names = TRUE)
	message("Top 5 predictions of the sample")
  for (i in 1:5){
    message(paste0("Predicted as ", names(agg_o[i]), " (" , round(agg_o[i]* 100, digits = 1), "%)"))
	}
  message(paste0("Prediction took ", round(as.numeric(difftime(time1 = b, time2 = a, units = "secs")), 3), " seconds"))
}
