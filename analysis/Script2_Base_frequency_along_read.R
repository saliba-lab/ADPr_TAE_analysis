###################### SCRIPT INFO ######################
# Script to quantify the occurence of each mutation (ATCG-) at every position in all reads.

# KNOWN LIMITATION: can't handle multiple samples with different reference sequences

# Set up directories ====
wd_dir <- "."
setwd(wd_dir)

# Load sample sheet (with check for correct function) 
sample_sheet <- read.csv2("Script2_Sample_sheet.csv", check.names = FALSE)

if (ncol(sample_sheet) == 1) { #To reload with read.csv if necessary
  sample_sheet <- read.csv("Script2_Sample_sheet.csv", check.names = FALSE)
}

# Setup output directory
output_folder <- "../outputs/Output_files_Script2"
dir.create(output_folder, recursive = TRUE)
setwd(output_folder)

# Locate input files
raw_data <- "../../data/Input_files_Script2"

# Run edition profiling ====
sample_files <- list.files(raw_data)[grep("\\.txt$", list.files(raw_data))]
edit_types <- c("A", "T", "C", "G", "-")

# Load all reference nucleotide sequences across every samples
all_ref_seq <- c() # initialize vector for all reference sequences
for (samp in sample_files) { #get through the different input files
  data <- read.delim(paste0(raw_data, "/", samp), check.names = FALSE)
  
  if (length(grep("-", data$Reference_Sequence)) == 0) {
    ref_seq <- unique(data$Reference_Sequence)
  } else {
    ref_seq <- data$Reference_Sequence[-grep("-", data$Reference_Sequence)]
    ref_seq <- unique(ref_seq)
  }
  
  all_ref_seq <- c(all_ref_seq, ref_seq)
  
} #end for loop samp

# Execute edition profiling
if ( length(unique(all_ref_seq)) == 1 ) { # Execute script only if all samples have the same reference sequence
  # Analysis
  table_list <- list() #initialize list to store all combined_result tables
  
  for (edit in edit_types) {
    combined_table <- data.frame() #initialize combined table dataframe
    sample_files <- list.files(raw_data)[grep("\\.txt$", list.files(raw_data))]
    
    for (samp in sample_files) { #get through the different input files
      data <- read.delim(paste0(raw_data, "/", samp), check.names = FALSE)
      thresh_counts <- sample_sheet$Threshold_read_counts[sample_sheet$Sample == samp]
      
      if (length(grep("-", data$Reference_Sequence)) == 0) {
        ref_seq <- unique(data$Reference_Sequence)
      } else {
        ref_seq <- data$Reference_Sequence[-grep("-", data$Reference_Sequence)]
        ref_seq <- unique(ref_seq)
      }
      
      combined_table <- rbind(combined_table, c("Reference_sequence", unlist(strsplit(ref_seq, ""))))
      edit_percent <- c() #initialize vector storing edition percentage
      
      # Loop through each position and count "-"
      for (base in 1:nchar(ref_seq)) {
        extracted_base <- substr(data$Aligned_Sequence, start = base, stop = base)
        indexes <- grep(edit, extracted_base)
        sub_data <- data[indexes, ]
        sub_data <- subset(sub_data, subset = sub_data$`#Reads` > thresh_counts)
        edit_percent <- c(edit_percent, round(sum(sub_data$`%Reads`), digits = 2))
      } #end for loop base
      
      combined_table <- rbind(combined_table, c(samp, edit_percent))
      
    } #end for loop samp
    
    combined_table[, 1] <- gsub(".txt", "", combined_table[, 1]) #remove .txt from the file names
    combined_table <- combined_table[!duplicated(combined_table), ] #remove multiple occurences of the reference sequence
    
    table_list[[edit]] <- combined_table # fill up the list with combined table for this edit type
  } # end for loop edit
  
  # Now calculate the total percentage of mutagenesis at each position in each sample
  percent_mut_table <- data.frame() # initialize table
  
  for (ref_base in c("A", "T", "C", "G")) {
    df <- table_list[[ref_base]]
    df <- apply(
      df[ ,-1],
      MARGIN = 2,
      function(x) {
        if (x[1] == ref_base) {
          x[2:length(x)] <- 0
        }
        return(x)
      }
    )
    df <- cbind(table_list[[ref_base]][, 1], df)
    percent_mut_table <- rbind(percent_mut_table, df)
  } #end for loop ref_base
  
  # Reshape table
  percent_mut_table <- percent_mut_table[percent_mut_table$V1 != "Reference_sequence", ]
  percent_mut_table[, -1] <- sapply(percent_mut_table[, -1], as.numeric, na.rm = TRUE) #convert to numerical values
  percent_mut_table <- aggregate(. ~ V1, data = percent_mut_table, sum, na.rm = TRUE) #sumup together the percentages for a given sample
  colnames(percent_mut_table) <- colnames(table_list[[1]])
  percent_mut_table <- rbind(table_list[[1]][1, ], percent_mut_table)
  table_list[["percent_mutagenesis"]] <- percent_mut_table #save the final table in the list
  
} else { # Abort and give warning if all samples do not share the same reference sequence
  print("ERROR: All samples do not have the same reference sequence")
}

# Save the table_list as an excel worksheet ====
openxlsx::write.xlsx(
  table_list, 
  "complete_MBE.xlsx", 
  rowNames = FALSE, 
  colNames = FALSE
)

# Save output table with replicates combined together ====
replicate_table_list <- list() #initialize table_list
edit_types <- c("A", "T", "C", "G", "-", "percent_mutagenesis")

# Extract group list and order by decreasing number of replicates, to start for loop with highest number of replicates
groups <- table(sample_sheet$Group)
group_list <- names(groups)[order(groups, decreasing = TRUE)]

for (samp_group in group_list) {
  group_tab <- data.frame( #initialize dataframe for sample group
    Reference_seq = unlist(strsplit(ref_seq, split = ""))
  ) 
  
  for (edit in edit_types) {
    #Select samples from the same replicate
    samples <- sample_sheet$Sample[sample_sheet$Group == samp_group]
    samples <- gsub(".txt", "", samples)
    
    
    tab <- table_list[[edit]][table_list[[edit]]$X.Reference_sequence. %in% samples, ]
    tab <- as.data.frame(t(tab))
    colnames(tab) <- paste0(edit, "_Rep", seq(from = 1, to = ncol(tab), by = 1))
    tab <- tab[-1, ]
    
    colnames(tab) <- gsub("^A", "Percent_A", colnames(tab))
    colnames(tab) <- gsub("^T", "Percent_T", colnames(tab))
    colnames(tab) <- gsub("^C", "Percent_C", colnames(tab))
    colnames(tab) <- gsub("^G", "Percent_G", colnames(tab))
    colnames(tab) <- gsub("^-", "Percent_deletion", colnames(tab))
    
    group_tab <- cbind(group_tab, tab)
  } #end for loop edit
  
  sample_idx <- paste(samples, collapse = ",")
  replicate_table_list[[sample_idx]] <- group_tab
  
} #end for loop samp_group

# Save the table_list as an excel worksheet ====
openxlsx::write.xlsx(
  replicate_table_list, 
  "complete_MBE_per_replicates.xlsx", 
  rowNames = FALSE, 
  colNames = TRUE
)

###################### END OF SCRIPT ######################