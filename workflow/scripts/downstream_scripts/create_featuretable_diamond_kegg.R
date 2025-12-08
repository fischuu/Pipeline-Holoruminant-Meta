#args <- commandArgs(trailingOnly = TRUE)
#project_folder <- args[1]
#result_folder <- args[2]

# Load the required library
library("data.table")

# Set the paths (set those for manual execution)
# project_folder <- "/scratch/project_2010176/metaG_groupAssembly/"
# result_folder <- "results/read_annotate/diamond/kegg"

# Get the available result files
result_files <- list.files(file.path(project_folder, result_folder), pattern="*.out")

# Initialize an empty list to store the results, a counts one per occasion, b counts weights (1/number of options)
all_counts_a <- list()
all_counts_b <- list()

# Now loop through all files
for (sample_run in 1:length(result_files)) {
  # Feedback for function
  cat("Started to summarize result file", result_files[sample_run], "at", date(), "...")
  
  # Set the current sample
  current_file <- result_files[sample_run]
  
  # Import the data
  data <- fread(file.path(project_folder, result_folder, current_file), header = FALSE, sep = "\t", stringsAsFactors = FALSE, showProgress = TRUE)
  
  # Filter the top hits only (first row for each read)
  if(max_target_seqs) {
    top_hits <- data[!duplicated(data$V1), ]
  } else {
    top_hits <- data
  }
  
  # Extract keggs (in case the leading part in front of : is needed, turn it off here...)
  keggs <- sub("^[^:]*:(.*)", "\\1", top_hits$V2)
  # keggs <- top_hits$V2
  
  # Create table for counts (table_counts_a)
  table_counts_a <- table(keggs)

  # Store in list
  all_counts_a[[current_file]] <- table_counts_a

  # Feedback of the function
  cat("done!\n")
}

# Create a list of all names across the vectors
all_names_a <- unique(unlist(lapply(all_counts_a, names)))

# Function to align each vector to the common set of names, filling with NA where necessary
align_to_names <- function(x, all_names) {
  # Create a vector of NAs for the common set of names
  result <- rep(NA, length(all_names))
  # Set the names of the result to the common set of names
  names(result) <- all_names
  # Assign values to the corresponding positions based on the names of the current vector
  result[names(x)] <- x
  return(result)
}

# Remove unnamed elements from each vector in all_counts_a
filtered_counts_a <- lapply(all_counts_a, function(x) {
  x <- x[nzchar(names(x))]  # Keep only elements with non-empty names
  return(x)
})

# Now align only the filtered vectors
aligned_counts_a <- lapply(filtered_counts_a, align_to_names, all_names = all_names_a)

cat("aligned_counts ready!\n")

# Combine the aligned vectors into a data frame
merged_counts_a <- rbind(aligned_counts_a[[1]], aligned_counts_a[[2]])
for(i in 3:length(aligned_counts_a)) merged_counts_a <- rbind(merged_counts_a, aligned_counts_a[[i]])

rownames(merged_counts_a) <- result_files

cat("rownames added ready!\n")

export_a <- t(merged_counts_a)

export_a <- data.frame(Feature = row.names(export_a), export_a)

write.table(export_a, file=file.path(project_folder, result_folder, "summary_a.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# We do not need the b output, but to fulfill the output requirements from snakeamke we writean identical file out
write.table(export_a, file=file.path(project_folder, result_folder, "summary_b.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# 
# # Now merge the pathway names into orthologs
# 
# procaryotes <- fread("/scratch/project_2010176/Tommi_Databases/KEGG/PROKARYOTES.DAT")
# 
# cat("procaryotes reading ready!\n")
# 
# # Ensure both are data.tables
# setDT(procaryotes)
# setDT(export_a)
# 
# # Merge KO into export based on Feature = gene
# export_merged <- merge(export_a, procaryotes, by.x = "Feature", by.y = "gene", all.x = TRUE)
# 
# cat("export_merged ready!\n")
# 
# # Remove Feature column (optional, if no longer needed)
#  export_merged[, Feature := NULL]
# 
# # Group by KO and sum all sample columns
# summary_by_KO <- export_merged[, lapply(.SD, sum, na.rm = TRUE), by = KO]
# 
# ko_info_1 <- fread("/scratch/project_2010176/Tommi_Databases/KEGG/KO00000")
# ko_info_2 <- fread("/scratch/project_2010176/Tommi_Databases/KEGG/KO00002")
# 
# cat("ko information read ready!\n")
# 
# export_merged_1 <- merge(ko_info_1, summary_by_KO, by.x = "KO", by.y = "KO", all = TRUE)
# export_merged_final <- merge(ko_info_2, export_merged_1, by.x = "KO", by.y = "KO", all = TRUE)
# 
# cat("last merging ready!\n")
# 
# export_merged_final[is.na(export_merged_final)] <- 0
# 
# # Export the data
# write.table(export_merged_final, file=file.path(project_folder, result_folder, "summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# #write.table(export_b, file=file.path(project_folder, result_folder, "summary_b.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# 
# cat("Complete!\n")
