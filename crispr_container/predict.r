#!/usr/bin/env Rscript
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3) # suppress TF messages
suppressPackageStartupMessages(suppressWarnings(library(deepG)))
library(deepG)
library(keras)
library(ggplot2)
library(optparse)
library(magrittr)
source("helpers.r")

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
    default = "/data/example.fasta"),
  make_option(c("-s", "--step"), type = "integer", default = 100,
help = "Step size to iterate though sequences [default %default]",
metavar = "number"),
  make_option(c("-g", "--gap"), type = "integer", default = 350,
help = "CRISPR gap size [default %default]",
metavar = "number"),
  make_option(c("-m", "--min_seq_len"), type = "integer", default = 150,
help = "min CRISPR size [default %default]",
metavar = "number"),
 make_option(c("-c", "--conf_cutoff"), type = "numeric", default = 0.8,
              help = "conf_cutoff [default %default]",
              metavar = "number"),
make_option(c("-p", "--pos_rate"), type = "numeric", default = 0.8,
              help = "pos_rate [default %default]",
              metavar = "number"),
  make_option(c("-b", "--batch_size"), type = "integer", default = 100,
help = "Number of samples processed in one batch [default %default]",
metavar = "number"));

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if output directory exists and create
if (!file.exists("/data/output")) {
  dir.create("/data/output")
}

a <- Sys.time()

model <- keras::load_model_hdf5("Ep.031-val_loss0.04-val_acc0.992.hdf5",
    compile = FALSE)

summary(model)
maxlen <- model$input_shape[[2]]

# changeme
fasta_index <- 1
return_summary <- TRUE
include_label <- FALSE

step <- 100 # opt
summary_method <- "mean"
crispr_gap <- 350 # opt
min_seq_len <- 150 # opt
conf_cutoff <- 0.8 # opt
pos_rate <- 0.8 # opt

conf_plot_path <- "/data/output/conf_plot.pdf"
crispr_list_path <- "/data/output/crispr_list.rds"
recalculate_margins <- FALSE

# load in files
fasta_df <- microseq::readFasta(opt$input)
nt_seq <- fasta_df$Sequence[fasta_index]

if (nchar(nt_seq) < maxlen) {
  stop(paste0("Input sequence must be at least ", maxlen, " characters long"))
}

pred <- bitmap_pred(fasta_df = NULL,
                    fasta_index = NULL,
                    batch_size = opt$batch_size,
                    return_summary = return_summary,
                    step = opt$step,
                    char_sequence = nt_seq,
                    model = model,
                    rc = FALSE,
                    range = NULL,
                    include_label = include_label,
                    summary_method = summary_method)

pred$conf_non_CRISPR <- 1 - pred$pred
names(pred) <- c("position", "conf_CRISPR", "conf_non_CRISPR")

p <- ggplot(pred, aes(x = position, y = conf_CRISPR)) + geom_point() + ylab("CRISPR confidence")

crispr_list <- filter_crispr(states_df = pred,
                             crispr_gap = opt$crispr_gap,
                             conf_cutoff = opt$conf_cutoff,
                             pos_rate = opt$pos_rate,
                             min_seq_len = opt$min_seq_len,
                             maxlen = maxlen)

if (is.null(crispr_list)) {
  print("No CRISPR candidates found")
} else {
  num_candidates <- length(crispr_list)
  start_index <- vector("numeric", num_candidates)
  end_index <- vector("numeric", num_candidates)
  for (i in 1:num_candidates) {
    start_index[i] <- min(crispr_list[[i]]$position)
    end_index[i] <- max(crispr_list[[i]]$position)
    
  }
  df <- data.frame(start_index, end_index)
}

if (recalculate_margins & length(crispr_list) > 0) {
  
  for (i in 1:nrow(df)) {
    
    new_start_pred <- df$start_index - 100
    pred_start <- bitmap_pred(batch_size = 1,
                              return_summary = FALSE,
                              step = 1,
                              char_sequence = nt_seq,
                              model = model,
                              range = c(new_start_pred, new_start_pred + maxlen - 1),
                              include_label = FALSE)
    
    new_end_pred <- df$end_index - 200
    pred_end <- bitmap_pred(batch_size = 1,
                            return_summary = FALSE,
                            step = 1,
                            char_sequence = nt_seq,
                            model = model,
                            range = c(new_end_pred, new_end_pred + maxlen - 1),
                            include_label = FALSE)
    
    df$start_index[i] <- pred_start$position[which(pred_start$pred > opt$conf_cutoff)] %>% min()
    df$end_index[i] <- pred_end$position[which(pred_end$pred > opt$conf_cutoff)] %>% max()
  }
}

if (!is.null(conf_plot_path)) ggsave(conf_plot_path, p)
if (!is.null(crispr_list_path)) saveRDS(df, crispr_list_path)