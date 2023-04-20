#!/usr/bin/env Rscript
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3) # suppress TF messages
suppressPackageStartupMessages(suppressWarnings(library(deepG)))
library(magrittr)
library(optparse)
library(ggplot2)
library(keras)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "/data/example.fasta"),
 make_option(c("-o", "--output"), type = "character", default = "/data/output.csv"),
  make_option(c("-b", "--batch_size"), type = "integer", default = 20,
help = "Number of samples processed in one batch [default %default]",
metavar = "number"),
 make_option(c("-s", "--step"), type = "integer", default = 1000,
help = "Step size to make prediction [default %default]",
metavar = "number"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

message("Load model")
model <- keras::load_model_hdf5("bacteria_spore_2023-01-23.hdf5", compile = FALSE)

message("Processing file")
pred <- predict_model(output_format = "one_seq",
                      model = model,
                      layer_name = "dense_1",
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

df <- data.frame(pred$states)
names(df) <- c("non-sporulating", "sporulating")

message("Writing predictions")
write.csv(df, file = opt$output, row.names = FALSE)

agg <- colMeans(df)
agg_o <- agg[order(agg, decreasing = T)]
agg_o <- agg_o[which(round(agg_o, digits = 2) > 0)]

for (i in 1:length(agg_o)){
    message(paste0(names(agg_o[i]), ": " , round(agg_o[i]* 100, digits = 1), "%"))
}
