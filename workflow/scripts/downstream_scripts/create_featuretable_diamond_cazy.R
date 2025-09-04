args <- commandArgs(trailingOnly = TRUE)
project_folder <- args[1]

# Load the required library
  library("data.table")

# Set the paths (set those for manual execution)
  # project_folder <- "/scratch/project_2010176/metaG_groupAssembly/"
  result_folder <- "results/preprocess/diamond/cazy"

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

# Create a list of all names across the vectors
  all_names_a <- unique(unlist(lapply(all_counts_a, names)))
  all_names_b <- unique(unlist(lapply(all_counts_b, names)))

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
  filtered_counts_b <- lapply(all_counts_b, function(x) {
    x <- x[nzchar(names(x))]  # Keep only elements with non-empty names
    return(x)
  })

# Now align only the filtered vectors
  aligned_counts_a <- lapply(filtered_counts_a, align_to_names, all_names = all_names_a)
  aligned_counts_b <- lapply(filtered_counts_b, align_to_names, all_names = all_names_b)


# Combine the aligned vectors into a data frame
  merged_counts_a <- rbind(aligned_counts_a[[1]], aligned_counts_a[[2]])
  merged_counts_b <- rbind(aligned_counts_b[[1]], aligned_counts_b[[2]])
  for(i in 3:length(aligned_counts_a)) merged_counts_a <- rbind(merged_counts_a, aligned_counts_a[[i]])
  for(i in 3:length(aligned_counts_b)) merged_counts_b <- rbind(merged_counts_b, aligned_counts_b[[i]])

  rownames(merged_counts_a) <- result_files
  rownames(merged_counts_b) <- result_files
  
#  export_a <- data.frame(sample_id=rownames(merged_counts_a), merged_counts_a)
#  export_b <- data.frame(sample_id=rownames(merged_counts_b), merged_counts_b)
#  colnames(export_a) <- c("sample_id", colnames(merged_counts_a))
#  colnames(export_b) <- c("sample_id", colnames(merged_counts_b))
#  
## Export the data
#  write.table(t(export_a), file=file.path(project_folder, result_folder, "summary_a.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
#  write.table(t(export_b), file=file.path(project_folder, result_folder, "summary_b.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  
  export_a <- t(merged_counts_a)
  export_b <- t(merged_counts_b)
  
  export_a <- data.frame(Feature = row.names(export_a), export_a)
  export_b <- data.frame(Feature = row.names(export_b), export_b)
  
  # Export the data
  write.table(export_a, file=file.path(project_folder, result_folder, "summary_a.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(export_b, file=file.path(project_folder, result_folder, "summary_b.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  
  cat("Complete!\n")