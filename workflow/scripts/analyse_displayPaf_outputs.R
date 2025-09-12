suppressMessages(library(dplyr))
suppressMessages(library(data.table))
library(ggplot2)
library("optparse")
library(stringr)
library(tidyr)

options(warn = -1)  # Suppress all warnings

# Define the option parser
parser <- OptionParser()
parser <- add_option(parser, c("-c", "--inputCsv"), type="character", help="displayPafAlignment CSV file")
parser <- add_option(parser, c("-o", "--outputPrefix"), type="character", default="<sample_id>_<region_id>", help="Output prefix for all output files [default %default]")

# Parse command-line arguments
args <- parse_args(parser)

# Print parsed arguments (for debugging)
print(args)

# Check if input file is provided
if (is.null(args$inputCsv)) {
  stop("Error: displayPafAlignment Csv file is required. Use -c or --inputCsv to specify the file", call. = FALSE)
}

# Define input directory
input_csv <- args$inputCsv
input_dir <- dirname(input_csv)

# Extract asm_preset from csv
  # Split by "_"
parts <- strsplit(basename(input_csv), split = "_")[[1]]
  # Extract "asmPreset" (2nd last element)
asm_version <- parts[length(parts) - 1]
  # Print the extracted value
print(asm_version)

# Define output directory
output_dir <- file.path(input_dir, paste0(asm_version,"_plots"))
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
# Define output prefix
output_prefix <- file.path(output_dir, args$outputPrefix)

# Define input files
displayPafCsv_file <- file.path(input_csv)

# Print paths for debugging
cat("Input Directory and input CSV is:", input_dir, displayPafCsv_file, "\n")
cat("Outputs will be dumped in:", output_dir, "\n")

df <- read.csv(displayPafCsv_file, sep = ",", col.names = c('Paf_line_number','Assembled','Reference','Alignment_position','Reference_position','Assembled_position','CIGAR'))

# Get unique references
unique_refs <- unique(df$Reference)

# Loop through each unique Reference and generate and save a plot
plots <- list()
for (ref in unique_refs) {
  df_subset <- df %>%
    filter(Reference == ref) %>%
    arrange(Reference_position) %>%
    mutate(row_number = row_number())  # Assign row number for the y-axis

  p <- ggplot(df_subset, aes(x = row_number, y = Reference_position)) +
    geom_text(aes(label = CIGAR, color = CIGAR), size = 3, fontface = "bold") +
    geom_line() +
    # geom_hline(yintercept = 1207458, linetype='dotted', col = '#000000') +  # For LGA sample
    # geom_hline(yintercept = 1240358, linetype='dotted', col = '#000000') +  # For LGA sample
    # geom_hline(yintercept = 1210823, linetype='dotted', col = '#000000') +  # For DEN63190 sample
    # geom_hline(yintercept = 1244807, linetype='dotted', col = '#000000') +  # For DEN63190 sample
    facet_grid(~Assembled) +
    labs(title = paste0("Assembly alignments to ", ref), x = "INDELs/Mismatches", y = "Reference Position") +
    theme_bw() +
    scale_color_manual(values = c("D" = "red", "I" = "blue", "X" = "#1d771d")) +
    theme(legend.position = "none")

  plots[[ref]] <- p

  # Save plot as PNG file
  filename <- paste0(output_prefix, "_varsPlot_", gsub("[:/]", "_", ref), ".png")  # Replace special characters in filenames
  ggsave(filename, plot = p, width = 12, height = 6, dpi = 300)
  paste0("Saves ", filename)
  paste0(" ")
}


# Loop through each unique Reference and generate and save a plot
plots <- list()
for (ref in unique_refs) {
  df_subset <- df %>%
    filter(Reference == ref) %>%
    arrange(Reference_position)
    
  p <- ggplot(df_subset, aes(x = Assembled_position, y = Reference_position)) +
    geom_point(aes(color = CIGAR), size = 1.5, alpha = 0.5) +
    geom_line() +
    # geom_hline(yintercept = 1207458, linetype='dotted', col = '#000000') +  # For LGA sample
    # geom_hline(yintercept = 1240358, linetype='dotted', col = '#000000') +  # For LGA sample
    # geom_hline(yintercept = 1210823, linetype='dotted', col = '#000000') +  # For DEN63190 sample
    # geom_hline(yintercept = 1244807, linetype='dotted', col = '#000000') +  # For DEN63190 sample
    facet_grid(~Assembled) +  # Facet by Assembled
    labs(title = paste0("Assembly alignments to ", ref), x = "Assembled Position", y = "Reference Position") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_color_manual(values = c("D" = "red", "I" = "blue", "X" = "#1d771d"))

  plots[[ref]] <- p

  # Save plot as PNG file
  filename <- paste0(output_prefix, "_plot_RvsA_", gsub("[:/]", "_", ref), ".png")  # Replace special characters in filenames
  ggsave(filename, plot = p, width = 12, height = 6, dpi = 300)
}

# Print plots
plots
