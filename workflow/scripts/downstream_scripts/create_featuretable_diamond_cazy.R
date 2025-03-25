# Load the required library
library("data.table")

# Set the paths
project_folder <- "/scratch/project_2001829/RuminomicsHighLow"
result_folder <- "results/preprocess/diamond/cazy"

# Get the available result files
result_files <- list.files(file.path(project_folder, result_folder), pattern="*.out")

# Initialize an empty list to store the results
all_counts_a <- list()
all_counts_b <- list()

# Now loop through all files
for (sample_run in 1:length(result_files)) {
  
#for (sample_run in 1:3) {
  
  # Feedback for function
  cat("Started to summarize result file", result_files[sample_run], "at", date(), "...")
  
  # Set the current sample
  current_file <- result_files[sample_run]
  
  # Import the data
  data <- fread(file.path(project_folder, result_folder, current_file), header = FALSE, sep = "\t", stringsAsFactors = FALSE, showProgress = FALSE)
  
  # Filter the top hits only (first row for each read)
  top_hits <- data[!duplicated(data$V1), ]
  
  # Extract CAZy families after the first pipe symbol
  cazy_families <- sub("^[^|]*\\|(.*)", "\\1", top_hits$V2)
  
  # Split entries by "|"
  split_labels <- unlist(strsplit(cazy_families, "\\|"))
  
  # Create table for counts (table_counts_a)
  table_counts_a_tmp <- table(split_labels)
  table_counts_a <- table_counts_a_tmp[!grepl(".fasta", names(table_counts_a_tmp))]
  
  # Remove ".fasta" entries from split list
  clean_fasta_entries <- function(lst) {
    lapply(lst, function(x) x[!grepl(".fasta", x)])
  }
  
  # Weighted count (table_counts_b)
  split_list <- strsplit(cazy_families, "\\|")
  cleaned_split_list <- clean_fasta_entries(split_list)
  weighted_counts <- unlist(lapply(cleaned_split_list, function(x) rep(1 / length(x), length(x)) * 1))
  names(weighted_counts) <- unlist(cleaned_split_list)
  table_counts_b <- tapply(weighted_counts, names(weighted_counts), sum)
  
  # Store in list
  all_counts_a[[current_file]] <- table_counts_a
  all_counts_b[[current_file]] <- table_counts_b
  
  # Feedback of the function
  cat("done!\n")
}

# Merge all count tables, filling missing values with NA
merged_counts_a <- as.data.frame(do.call(cbind, lapply(all_counts_a, function(x) setNames(as.numeric(x), names(x)))))
merged_counts_b <- as.data.frame(do.call(cbind, lapply(all_counts_b, function(x) setNames(as.numeric(x), names(x)))))

# Add row names
rownames(merged_counts_a) <- unique(unlist(lapply(all_counts_a, names)))
rownames(merged_counts_b) <- unique(unlist(lapply(all_counts_b, names)))

# Rename columns with sample names
colnames(merged_counts_a) <- names(all_counts_a)
colnames(merged_counts_b) <- names(all_counts_b)

# Compute within-sample correlation
within_sample_correlation <- sapply(1:length(result_files), function(i) {
  common_labels <- intersect(rownames(merged_counts_a), rownames(merged_counts_b))
  cor(merged_counts_a[common_labels, i], merged_counts_b[common_labels, i], use = "pairwise.complete.obs")
})

# Combine results into a single data frame
final_results <- list(
  raw_counts = merged_counts_a,
  weighted_counts = merged_counts_b,
  correlation = within_sample_correlation
)
