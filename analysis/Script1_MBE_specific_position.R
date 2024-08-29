###################### SCRIPT INFO ######################
# Script to quantify the total percentage of reads with substitution or deletion at a specific position
# Return 2x tables: 
# (1) First one with each row corresponding to one sample -> Combined_table.csv/xlsx
# (2) a second one where replicates of the same condition are grouped together on a single row -> Replicate_table.csv/xlsx
#
# Sample information regarding the replicates (group) and parameters must be reported in the sample sheet.
#
# KNOWN LIMITATION: can't handle multiple samples with different reference nucleotides

# Set up directories ====
wd_dir <- "."
setwd(wd_dir)

sample_sheet <- read.csv("Script1_Sample_sheet.csv", check.names = FALSE)

output_folder <- "../outputs/Output_files_Script1"
dir.create(output_folder, recursive = TRUE)
setwd(output_folder)

raw_data <- "../../data/Input_files_Script1"

# Run script ====
combined_table <- data.frame() #initialize combined result table
sample_files <- list.files(raw_data)[grep("\\.txt$", list.files(raw_data))]

for (samp in sample_files) {
  
  # Load data and parameters for selected sample
  data <- read.delim(paste0(raw_data, "/", samp), check.names = FALSE)
  mutation_pos <- sample_sheet$Mutation_position[sample_sheet$Sample == samp]
  thresh_counts <- sample_sheet$Threshold_read_counts[sample_sheet$Sample == samp]
  
  # Remove reference sequence with deletions
  if (length(grep("-", data$Reference_Sequence)) == 0) {
    ref_seq <- unique(data$Reference_Sequence)
  } else {
    ref_seq <- data$Reference_Sequence[-grep("-", data$Reference_Sequence)]
    ref_seq <- unique(ref_seq)
  }
  
  # Set reference sequence and decompose it per nucleotide
  ref_nt <- substr(ref_seq, start = mutation_pos, stop = mutation_pos)
  
  # Perform mutational analysis
  data$MBE <- substr(data$Aligned_Sequence, start = mutation_pos, stop = mutation_pos) # extract nucleotide in target position
  data <- subset(data, subset = MBE != ref_nt & data$`#Reads` > thresh_counts) # remove unedited read
  
  possible_MBE <- c("A", "T", "C", "G", "-") # list all possible nucleotides
  possible_MBE <- possible_MBE[-grep(ref_nt, possible_MBE)] # remove reference nucleotide from the list
  
  # Create output table
  results <- data.frame(
    Sample_name = gsub(".txt", "", samp),
    Percent_mut1 = sum(subset(data, subset = MBE == possible_MBE[1])$`%Reads`),
    Percent_mut2 = sum(subset(data, subset = MBE == possible_MBE[2])$`%Reads`),
    Percent_mut3 = sum(subset(data, subset = MBE == possible_MBE[3])$`%Reads`),
    Percent_total_substitution = sum(subset(data, subset = MBE != "-")$`%Reads`), # total is excluding deletions ??
    Percent_deletion = sum(subset(data, subset = MBE == "-")$`%Reads`)
  )
  
  column_names <- c( #define column names for the result table
    paste0("Percent_", possible_MBE[1]),
    paste0("Percent_", possible_MBE[2]),
    paste0("Percent_", possible_MBE[3])
  )
  
  colnames(results)[2:4] <- column_names #update column names
  
  combined_table <- rbind(combined_table, results)
  
} # end for loop


# Save output table with each sample separated ====
write.csv(combined_table, "Combined_table.csv", row.names = F) #save as csv file
openxlsx::write.xlsx(combined_table, "Combined_table.xlsx", rowNames = FALSE) #save as xlsx file for GraphPad

# Save output table with replicates combined together ====
replicate_table <- data.frame() #initialize table

# Extract group list and order by decreasing number of replicates, to start for loop with highest number of replicates
groups <- table(sample_sheet$Group)
group_list <- names(groups)[order(groups, decreasing = TRUE)]

for (samp_group in group_list) {
  #Select samples from the same replicate
  samples <- sample_sheet$Sample[sample_sheet$Group == samp_group]
  samples <- gsub(".txt", "", samples)
  
  
  tab <- combined_table[combined_table$Sample_name %in% samples, ]
  tab <- as.data.frame(
    c(
      data.frame(
        Samples = paste(tab$Sample_name, collapse = ", ")
      ),
      unlist(tab[, -1])
    )
  )
  
  #Handle case with only 1 replicate
  if (ncol(tab) == 6) {
    colnames(tab)[-1] <- paste0(
      colnames(tab)[-1],
      "1"
    )
  }
  
  
  replicate_table <- dplyr::bind_rows(replicate_table, tab)
} #end for loop

# Save Replicate tables
write.csv(replicate_table, "Replicate_table.csv", row.names = F) #save as csv file
openxlsx::write.xlsx(replicate_table, "Replicate_table.xlsx", rowNames = FALSE, keepNA = TRUE) #save as xlsx file for GraphPad

###################### END OF SCRIPT ######################
