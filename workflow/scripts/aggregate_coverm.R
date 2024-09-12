#!/usr/bin/env Rscript

library(tidyverse)
library(argparse)

parser <- ArgumentParser()

parser$add_argument(
  "-i", "--input-folder",
  type = "character",
  dest = "input_folder",
  help = "Folder containing the *.tsv files"
)

parser$add_argument(
  "-o", "--output-file",
  type = "character",
  dest = "output_file",
  help = "Output TSV file"
)

args <- parser$parse_args()
input_folder <- args$input_folder
output_file <- args$output_file
output_folder <- dirname(output_file)

dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

files <- list.files(args$input_folder, pattern = "*.tsv", full.names = TRUE)

sample_names <-
  files %>%
  basename() %>%
  str_remove(".tsv")

files %>%
  set_names(sample_names) %>%
  map(
    function(x) {
      read_tsv(
        file = x, col_types = cols(), col_names = c("sequence_id", "counts"),
        skip = 1
      )
    }
  ) %>%
  bind_rows(.id = "sample_id") %>%
  pivot_wider(
    names_from = "sample_id", values_from = "counts", values_fill = NA
  ) %>%
  write_tsv(output_file)
