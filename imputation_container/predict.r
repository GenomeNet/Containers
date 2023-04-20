#!/usr/bin/env Rscript
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3) # suppress TF messages
suppressPackageStartupMessages(suppressWarnings(library(deepG)))
library(magrittr)
library(optparse)
library(ggplot2)
library(keras)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "/data/example.fasta"),
  make_option(c("-o", "--output"), type = "character", default = "/data/output.fasta"),
  make_option(c("-t", "--threshold"), type = "integer", default = 0.5,
  help = "Threshold [default %default]", metavar = "number"),
  make_option(c("-b", "--batch_size"), type = "integer", default = 2,
help = "Number of samples processed in one batch [default %default]",
metavar = "number"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if output directory exists and create
if (!file.exists("/data/output")) {
  dir.create("/data/output")
}

args <- commandArgs(trailingOnly = TRUE)
message("Load model")

model <- keras::load_model_hdf5("model_imputation_maxlen100.hdf5", compile = FALSE)

sequence <- microseq::readFasta(opt$input)$Sequence
len <- nchar(sequence)

maxlen <- 100
message("Processing FASTA file")
pred <- predict_model(vocabulary = c("a", "c", "g", "t", "n"),
                      output_format = "one_seq",
                      model = model,
                      layer_name = "flatten_4",
                      sequence = sequence,
                      round_digits = 4,
                      filename = NULL,
                      step = maxlen,
                      batch_size = 512,
                      verbose = FALSE,
                      return_states = TRUE,
                      padding = "standard",
                      mode = "label",
                      format = "fasta")


if (len > maxlen) {
  pred1 <- predict_model(vocabulary = c("a", "c", "g", "t", "n"),
                      output_format = "one_seq",
                      model = model,
                      layer_name = "flatten_4",
                      sequence = substr(sequence, len - maxlen + 1, len),
                      round_digits = 4,
                      filename = NULL,
                      step = maxlen,
                      batch_size = 4,
                      verbose = FALSE,
                      return_states = TRUE,
                      padding = "standard",
                      mode = "label",
                      format = "fasta")
} else{
  k <- (maxlen - len) * 4
}

num_of_pred <- len %/% maxlen
last_pred <- num_of_pred * maxlen
last_piece <- len %% maxlen
last_addition <- (maxlen - last_piece) * 4
n_positions <- which(strsplit(sequence, "")[[1]] == "N")
orig_ambigous_n <- length(which(strsplit(sequence, "")[[1]] == "N"))
message(paste0("Found ", length(n_positions), " ambigous nucleotides in the file"))

theshold <- opt$threshold

i <- 1
if (len > maxlen) {
  while (i <= length(n_positions)) {
    j <- n_positions[i]
    j_rem_maxlen <- ((j - 1) %% maxlen) + 1
    if  (j <= last_pred) {
        a <- which.max(pred$state[((j - 1) %/% maxlen) + 1,(j_rem_maxlen * 4 - 3):(j_rem_maxlen * 4)])
        a_prob <- max(pred$state[((j - 1) %/% maxlen) + 1,(j_rem_maxlen * 4 - 3):(j_rem_maxlen * 4)])

      if (a == 1){
          
        if (a_prob >= theshold) {
            substring(sequence, j) <- "A"
             message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'A' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
        } else {
                 message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
        }
          
      } else if (a == 2){
      
        if (a_prob >= theshold) {
            substring(sequence, j) <- "C"
                            message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'C' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
        } else {
                 message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
        }
          
      } else if (a == 3){
       
        if (a_prob >= theshold) {
            substring(sequence, j) <- "G"
                            message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'G' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
        } else {
                 message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
        }
          
      } else if (a == 4){
        
        if (a_prob >= theshold) {
            substring(sequence, j) <- "T"
                            message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'A' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
        } else {
                 message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
        }
          
      } 
    } else {
      a = which.max(pred1$state[1,(last_addition + 
                                   j_rem_maxlen * 4 - 3):(last_addition + j_rem_maxlen * 4)])
      a_prob = which.max(pred1$state[1,(last_addition + 
                                        j_rem_maxlen * 4 - 3):(last_addition + j_rem_maxlen * 4)])

      if (a == 1){
           
            if (a_prob >= theshold) {
                substring(sequence, j) <- "A"
                  message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'A' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
            } else {
                 message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
            }
          
          
      }else if (a == 2){
       
          
            if (a_prob >= theshold) {
                substring(sequence, j) <- "C"
                   message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'C' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
            } else {
                 message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
            }

      }else if (a == 3){
          
          if (a_prob >= theshold) {
            substring(sequence, j) <- "G"
               message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'G' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
        } else {
                 message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
        }
          
      }else if (a == 4){

          if (a_prob >= theshold) {
            substring(sequence, j) <- "T"
               message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'T' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
        } else {
                 message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
        }
      }
    } 
  i <- i + 1
  }
} else {
  while (i <= length(n_positions)) {
    j <- n_positions[i]
    j_rem_maxlen <- ((j - 1) %% maxlen) + 1
    a <- which.max(pred$state[1,(k + j_rem_maxlen * 4 - 3):(k + j_rem_maxlen * 4)])
    a_prob <- which.max(pred$state[1,(k+j_rem_maxlen * 4 - 3):(k+j_rem_maxlen * 4)])
    if (a == 1){
         if (a_prob >= theshold) {
               substring(sequence, j) <- "A"
                 message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'A' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
         } else {
           message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
         }
    
    }else if (a == 2) {
            if (a_prob >= theshold) {
               substring(sequence, j) <- "C"
                 message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'C' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
         } else {
           message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
         }
        
    }else if (a == 3) {
           if (a_prob >= theshold) {
               substring(sequence, j) <- "G"
                 message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'G' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
         } else {
           message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
         }
        
    }else if (a == 4) {
            if (a_prob >= theshold) {
               substring(sequence, j) <- "T"
                 message(paste0("Position " , j, ": Imputing nucleotide 'N' with 'T' (probability ",
                       round(a_prob, digits = 2) * 100, "%)"))
         } else {
           message(paste0("Position " , j, ": Skipping ambigous nucleotide 'N' since no probability is above threshold"))
         }
        
    }
  i <- i + 1
  }
}

seq <- microseq::readFasta(opt$input)
seq$Sequence <- sequence
microseq::writeFasta(seq, opt$output)

new_ambigous_n <- length(which(strsplit(sequence, "")[[1]] == "N"))

message(paste0("Successfully imputed ", orig_ambigous_n - new_ambigous_n,
               " positions out of ", orig_ambigous_n, " that meet criteria"))