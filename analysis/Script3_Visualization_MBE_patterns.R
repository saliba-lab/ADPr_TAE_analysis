###################### SCRIPT INFO ######################
# Script to quantify reads with multiple mutations and create a graphical representation of the sequences.
# Current script can handle mutations at 2 different positions as specified by the user in the sample sheet.
# The script generates matched sequence plots and bar plots; the bar plots indicate the abundancy of each sequence represented in the sequence plots.

# Set up directories ====
wd_dir <- "."
setwd(wd_dir)

# Load sample sheet (with check for correct function) 
sample_sheet <- read.csv2("Script3_Sample_sheet.csv", check.names = FALSE)

if (ncol(sample_sheet) == 1) { #To reload with read.csv if necessary
  sample_sheet <- read.csv("Script3_Sample_sheet.csv", check.names = FALSE)
}

backup_sample_sheet <- sample_sheet

# Set up output directory
output_folder <- "../outputs/Output_files_Script3"
dir.create(output_folder, recursive = TRUE)
setwd(output_folder)

# Setup plot directory
graph_dir <- "Plots"
dir.create(graph_dir)

# Locate input files
raw_data <- "../../data/Input_files_Script3"

# Get list of sample files ====
sample_files <- list.files(raw_data)[grep("\\.txt$", list.files(raw_data))]

if (sum(!sample_files %in% sample_sheet$Sample) > 0) {
  print("INFO: Some samples are missing from the sample sheet and won't be considered for analysis.")
  write.csv(
    sample_files[!sample_files %in% sample_sheet$Sample], 
    paste0("INFO_Samples_absent_from_samplesheet.csv")
  )
  
  sample_files <- sample_files[sample_files %in% sample_sheet$Sample]
}

if (sum(!sample_sheet$Sample %in% sample_files) > 0) {
  print("INFO: Some samples indicated in the sample sheet are not present in Input_files.")
  write.csv(
    sample_sheet$Sample[!sample_sheet$Sample %in% sample_files],
    paste0("INFO_Samples_absent_from_Input_files.csv")
  )
  
  sample_sheet <- sample_sheet[sample_sheet$Sample %in% sample_files, ]
}

# Iterate through each sample file ====
for (samp in sample_files) {
  
  # Read the data from the file and parameters
  data <- read.delim(paste0(raw_data, "/", samp), check.names = FALSE)
  data <- data[-grep("-", data$Reference_Sequence), ]
  
  mutation_positions <- c(
    sample_sheet$Mutation_position1[sample_sheet$Sample == samp],
    sample_sheet$Mutation_position2[sample_sheet$Sample == samp]
  )
  
  thresh_counts <- sample_sheet$Threshold_read_counts[sample_sheet$Sample == samp]
  
  # Initialize temporary and final selected rows
  final_selected_rows <- data.frame()
  
  # Select rows based on the first mutation position
  first_pos <- mutation_positions[1]
  ref_nt1 <- substr(data$Reference_Sequence, start = first_pos, stop = first_pos)
  data$MBE1 <- substr(data$Aligned_Sequence, start = first_pos, stop = first_pos)
  
  # Filter rows where the character at the first mutation position is different
  selected_rows_temp <- subset(data, MBE1 != ref_nt1 & data$`#Reads` > thresh_counts)
  
  # Iterate through the remaining mutation positions
  for (pos in mutation_positions[-1]) {
    ref_nt <- substr(selected_rows_temp$Reference_Sequence, start = pos, stop = pos)
    selected_rows_temp$MBE <- substr(selected_rows_temp$Aligned_Sequence, start = pos, stop = pos)
    
    # Filter rows where the character at the current mutation position is different
    selected_rows <- subset(selected_rows_temp, MBE != ref_nt)
    
    # Append the selected rows to the final result
    final_selected_rows <- rbind(final_selected_rows, selected_rows)
  }
  
  # Check if there are any rows in the final result
  if (nrow(final_selected_rows) == 0) {
    final_selected_rows <- data.frame(Message = "No data found matching the criteria")
  }
  
  # Write the final selected rows to an Excel file with the same name as the input file
  output_file <- paste0(sub(".txt", "", samp), "_output.xlsx")
  openxlsx::write.xlsx(final_selected_rows, output_file, rowNames = FALSE)
} # end samp for loop

# Visualization - PER REPLICATE ====
# Reload sample_sheet
sample_sheet <- backup_sample_sheet

# Load data and assemble tables per sample group
group_list <- list() #initialize table list
output_files <- list.files(getwd())[grep("\\.xlsx$", list.files(getwd()))]

# Filter out empty output tables
empty_files <- sapply( 
  output_files, function(tab) {
    tab <- openxlsx::read.xlsx(tab)
    tab <- tab[1, 1] == "No data found matching the criteria"
  }
)

if (sum(empty_files) > 0) {
  print("INFO: Some samples do not have mutated reads matching criteria and were excluded from further analysis.")
  write.csv(
    output_files[empty_files], 
    paste0("INFO_No_matching_criteria_samples.csv")
  )
}

output_files <- output_files[!empty_files]
sample_sheet <- sample_sheet[sample_sheet$Sample %in% gsub("_output.xlsx", ".txt", output_files), ]


# Assemble tables per sample group
for (current_group in unique(sample_sheet$Group)) { 
  group_samples <- sample_sheet$Sample[sample_sheet$Group == current_group]
  group_samples <- gsub(".txt", "_output.xlsx", group_samples)
  group_table <- data.frame()#intialize group table
  
  for (samp in group_samples) {
    
    tab <- openxlsx::read.xlsx(paste0(samp))
    tab$Group <- current_group
    tab$Sample <- samp
    tab$Sample <- gsub("_output.xlsx", "", tab$Sample)
    
    group_table <- rbind(group_table, tab)
  }
  
  group_list[[current_group]] <- group_table
} # end current_group for loop

# Filter deletions and add Seq_id
group_list <- lapply( # Remove reads with deletions at mutation positions
  group_list, function(tab) {
    
    index_no_deletion <- intersect(
      which(tab$MBE1 != "-"),
      which(tab$MBE != "-")
    )
    
    tab <- tab[index_no_deletion, ]
  }
) # end lapply

# Filter out sample groups having all mutated reads involving deletions
if (sum(sapply(group_list, function(x) dim(x)[1] == 0)) > 0) {
  print("INFO: Some groups of samples had only mutated reads involving deletions. These groups of samples were excluded from further analysis.")
  group_list <- group_list[sapply(group_list, function(x) dim(x)[1]) > 0]
  
  write.csv(
    unique(sample_sheet$Group)[!unique(sample_sheet$Group) %in% names(group_list)] , 
    paste0("INFO_Sample_groups_with_only_deletions.csv")
  )
}

group_list <- lapply( # Add Seq_id. Identical Aligned sequences will have the same Seq_id across replicates.
  group_list, function(tab) {
    
    # Create Seq_id lookup table
    lookup_Seq_id <- data.frame(
      Aligned_seq = unique(tab$Aligned_Sequence),
      Seq_id = paste0("Seq_", 1:length(unique(tab$Aligned_Sequence)))
    )
    
    # Add Seq_id
    tab$Seq_id <- sapply(
      tab$Aligned_Sequence, function(id) {
        lookup_Seq_id$Seq_id[match(id, lookup_Seq_id$Aligned_seq)]
      }
    )

    tab
  }
)

# Reshape tables for nucleotide sequence plotting
group_nt_list <- group_list

# Filter duplicated sequences
group_nt_list <- lapply(
  group_nt_list, function(tab) {
    tab <- tab[!duplicated(tab$Seq_id), ]
    tab
  }
)

# Reshape
group_nt_list <- lapply(
  group_nt_list, function(tab) {
    
    # Extract reference sequence
    ref_seq <- unlist(strsplit(unique(tab$Reference_Sequence), split = ""))
    
    # Initialize table with reference seq
    data <- data.frame( 
      Seq_id = "Reference_seq",
      Position = c(1:length(ref_seq)),
      Base = ref_seq
    )
    
    # Insert blank row for plotting
    data <- rbind(
      data,
      data.frame( 
        Seq_id = "",
        Position = c(1:length(ref_seq)),
        Base = "empty"
      )
    )
    
    # Reshape table
    for (seq_row in 1:nrow(tab)) {
      
      current_seq <- tab[seq_row, ]
      aligned_seq <- unlist(strsplit(current_seq$Aligned_Sequence, split = ""))
      
      tmp <- data.frame(
        Seq_id = rep(current_seq$Seq_id, length(ref_seq)),
        Position = c(1:length(ref_seq)),
        Base = aligned_seq
      )
      
      data <- rbind(data, tmp)
    }# end seq_row for loop
    data
  }
)

# Plotting nucleotide sequences
nucleotide_palette <- c(
  "A" = "#009E73",
  "C" = "#0072B2",
  "G" = "#000000",
  "T" = "#D55E00",
  "-" = "#999999",
  "empty" = "white"
)

sequence_plots <- lapply(
  group_nt_list, function(data) {
    
    data$Seq_id <- factor(data$Seq_id, levels = rev(unique(data$Seq_id)))
    
    ggplot2::ggplot(data, ggplot2::aes(x = Position, y = Seq_id, fill = Base)) +
      ggplot2::geom_tile(
        color = "white",
        lwd = 1.5,
        linetype = 1
      ) +
      ggplot2::scale_fill_manual(values = nucleotide_palette) +
      ggplot2::geom_text(
        ggplot2::aes(label = Base), 
        color = "white", 
        size = 4
      ) +
      ggplot2::coord_fixed() +
      cowplot::theme_cowplot(12) +
      ggplot2::theme(
        legend.position = "none",
        axis.title.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank()
      )
    
  }
)

# Save sequence plots
dir.create("Plots/Sequence_plots", recursive = TRUE)
graph_dir <- "Plots/Sequence_plots/"

for (group in names(sequence_plots)) {
  
  samples <- sample_sheet$Sample[sample_sheet$Group == group]
  samples <- gsub(".txt", "", samples)
  samples <- paste0(samples, collapse = "__")
  
  filename <- paste0(graph_dir, samples, "_sequence_plot.eps")
  ggplot2::ggsave(
    filename = filename,
    plot     = sequence_plots[[group]],
    width    = 14,
    height   = 7
  )
  
  filename <- paste0(graph_dir, samples, "_sequence_plot.pdf")
  ggplot2::ggsave(
    filename = filename,
    plot     = sequence_plots[[group]],
    width    = 14,
    height   = 7
  )
  
} # end group for loop

# Plotting read percentage per sequence
bar_plots <- lapply(
  group_list, function(data) {
    
    data$Seq_id <- factor(data$Seq_id, levels = rev(unique(data$Seq_id)))
    
    ggplot2::ggplot(data, ggplot2::aes(x = `%Reads`, y = Seq_id)) +
      #ggplot2::geom_col(stat = "count") +
      ggplot2::stat_summary(geom = "col", fun = mean, fill = "lightgrey") +
      ggplot2::geom_point(
        data = data, 
        ggplot2::aes(
          x = `%Reads`, 
          y = Seq_id, 
          col = Sample
        )
      ) +
      ggplot2::ylab("") +
      ggplot2::xlab("Percentage of total reads (%)") +
      ggplot2::scale_x_continuous(
        expand = c(0, 0), 
        limits = c(0, max(data$`%Reads`) + 0.05), 
        position = "top"
      ) +
      cowplot::theme_cowplot(12)
  }
)

# Save bar plots
dir.create("Plots/Bar_plots", recursive = TRUE)
graph_dir <- "Plots/Bar_plots/"

for (group in names(sequence_plots)) {
  
  samples <- sample_sheet$Sample[sample_sheet$Group == group]
  samples <- gsub(".txt", "", samples)
  samples <- paste0(samples, collapse = "__")
  
  filename <- paste0(graph_dir, samples, "_bar_plot.eps")
  ggplot2::ggsave(
    filename = filename,
    plot     = bar_plots[[group]],
    width    = 14,
    height   = 7
  )
  
  filename <- paste0(graph_dir, samples, "_bar_plot.pdf")
  ggplot2::ggsave(
    filename = filename,
    plot     = bar_plots[[group]],
    width    = 14,
    height   = 7
  )
  
}
###################### END OF SCRIPT ######################