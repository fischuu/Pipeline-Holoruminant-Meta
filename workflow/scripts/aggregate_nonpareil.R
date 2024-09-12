#!/usr/bin/Rscript

library(tidyverse)
library(argparse)
library(Nonpareil)

parser <- ArgumentParser()

parser$add_argument(
  "-i", "--input-folder",
  type = "character",
  dest = "input_folder",
  help = "Folder containing the *.npo files"
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

nonempty_files <-
  list.files(args$input_folder, pattern = "*.npo", full.names = TRUE) %>%
  data.frame(file = .) %>%
  mutate(size = file.info(file)$size) %>%
  filter(size > 0) %>%
  pull(file)

nonempty_files %>%
  Nonpareil.set(plot = FALSE) %>%
  summary() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  as_tibble() %>%
  write_tsv(output_file)
