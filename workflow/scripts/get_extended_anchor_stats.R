suppressMessages(library(dplyr))
suppressMessages(library(data.table))
library(ggplot2)
library("optparse")
library(stringr)
library(tidyr)
library(patchwork)

options(warn = -1)  # Suppress all warnings

# Define the option parser
parser <- OptionParser()
# parser <- add_option(parser, c("-a", "--anchorsJsonTsv"), type="character", help="Anchors to read info TSV (output of preprocess_vganchor_outfiles.py)")
# parser <- add_option(parser, c("-b", "--anchorsInfoTsv"), type="character", help="Anchor sizes TSV. It contains information on all anchors found from the pangenome.")
parser <- add_option(parser, c("-d", "--inputDir"), type="character", help="Directory with anchor results for the region of interest")
parser <- add_option(parser, c("-o", "--outputPrefix"), type="character", default="<sample_id>_<region_id>", help="Output prefix for all output files [default %default]")
parser <- add_option(parser, c("-r", "--region"), type="character", help="Provide coordinates of the subregion in format chrom:start-end")
# parser <- add_option(parser, c("-t", "--tech"), type="character", help="Sequencing technology")


# Parse command-line arguments
args <- parse_args(parser)

# Print parsed arguments (for debugging)
print(args)

# Check if input file is provided
if (is.null(args$inputDir)) {
  stop("Error: Anchors directory name is required. Use -d or --inputDir to specify the directory.", call. = FALSE)
}

# Define input directory
input_dir <- args$inputDir
# Define output directory
output_dir <- file.path(input_dir, "extended_anchor_stats")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
# Define output prefix
output_prefix <- file.path(output_dir, args$outputPrefix)
# parse region
matches <- str_match(args$region, "(.*):(.*)-(.*)")
chr <- matches[2]
start <- as.integer(matches[3])
end <- as.integer(matches[4])
cat("region of interest is:", "\n")
cat("Chromosome:", chr, "\n")
cat("Start:", start, "\n")
cat("End:", end, "\n")


# Define input files
anchorsJsonTsv_file <- file.path(input_dir, "anchors/extended_anchor_reads_info.tsv")
anchorsInfoTsv_file <- file.path(input_dir, "anchors/subgraph.anchors.json.subgraph.sizes.extended.tsv")
shastaAnchorsJson_file <- file.path(input_dir, "shasta/ShastaRun/AnchorsFromJson.csv")
subgraphGaf_file <- file.path(input_dir, "chunk/subgraph.gaf")
read_processed_tsv_file <- file.path(input_dir, "anchors/subgraph.anchors.json.reads_processed.tsv")

# Print paths for debugging
cat("Input Directory is:", input_dir, "\n")
cat("Outputs will be dumped in:", output_dir, "\n")

anchor_read_stats <- read.csv(anchorsJsonTsv_file, sep = "\t", col.names = c('Anchor_path', 'Anchor_length', 'bp_matched_reads'))
anchor_sizes <- read.csv(anchorsInfoTsv_file, sep = "\t", header = TRUE)
anchor_sizes <- anchor_sizes %>% group_by(Anchor_path) %>% mutate(nodes_matching_reads=length(unlist(strsplit(Paths_associated_with_anchor, ","))))
anchor_shasta_id <- read.csv(shastaAnchorsJson_file, sep = ",", col.names = c('Anchor_path', 'forward_id', 'reverse_id'))

paste0("Total number of anchors found in the pangenome for region ", args$region, " are ", nrow(anchor_sizes))
paste0("Total number of anchors provided to shasta (i.e. anchors with one or more bp matched reads) for region ", args$region, " are ", nrow(anchor_read_stats))

# Outer merge both dataframes
anchor_stats <- merge(x=anchor_sizes, y=anchor_read_stats, by = c("Anchor_path", "Anchor_length", "bp_matched_reads"), all=TRUE)
anchor_stats$bp_matched_reads[is.na(anchor_stats$bp_matched_reads)] <- 0

# merge with Shasta's AnchorsFromJson.csv
anchor_stats_all <- merge(x=anchor_stats, y=anchor_shasta_id, by = "Anchor_path", all=TRUE)

anchor_stats_all <- anchor_stats_all %>%
        mutate(snarl_id=as.character(snarl_id)) %>%
        mutate(chrom=chr,
               start=Anchor_pos_in_ref_path,
               end=Anchor_pos_in_ref_path+Anchor_length)

## Remove anchors with a 0 Anchor_pos_in_ref_path and output that in the log
## TODO: Handle this better. Estimate Anchor_pos if that anchor isn't present in CHM13.
anchor_stats_with_ref_paths_only <- anchor_stats_all %>% filter(!is.na(Anchor_pos_in_ref_path) & Anchor_pos_in_ref_path != 0)
anchor_stats_with_ref_paths_only <- anchor_stats_with_ref_paths_only %>% mutate(snarl_id=as.character(snarl_id))
paste0("Total number of anchors given to shasta for region (after removing Anchor_pos_in_ref_path=0 anchors)", args$region, " are ", nrow(anchor_stats_with_ref_paths_only))


###################### GET STATS AND PLOTS ############################
####### 1. Anchor length distribution
lens <- anchor_stats_all %>% 
  ggplot(aes(x = Anchor_length)) + 
  geom_histogram(binwidth = 10, fill = "#cddd9c", color="black") +
  geom_vline(xintercept = 20, linetype = "longdash", colour = "red") +
  # geom_vline(xintercept = 40, linetype = "longdash", colour = "blue") +
  theme_bw() + 
  labs(x = "Anchor length", y = "Number of Anchors", title = "Distribution of anchor lengths") +
  scale_x_continuous(breaks =seq(0,max(anchor_stats_all$Anchor_length)+100,100)) +
#   scale_y_continuous(limits = c(0,120)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 0, vjust = 0.5))

f1 <- paste0(output_prefix, "_anchor-lengths.png")
ggsave(filename=f1,plot=lens,width=10,height=6,dpi=300)
paste("")
paste0("Generated ", f1)


###### 2. Read coverage of anchors distribution
# bp_matched_reads distribution
read_c <- anchor_stats_all %>% 
  ggplot(aes(x = bp_matched_reads)) + 
  geom_histogram(binwidth = 1, fill = "#9cc9dd", color="black") +
  theme_bw() + 
  labs(x = "Read coverage", y = "Number of Anchors", title = "Number of Anchors per read depth") +
#   scale_x_continuous(breaks =seq(0,100,1)) +
#   scale_y_continuous(limits = c(0,120)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

f2 <- paste0(output_prefix, "_read-counts.png")
ggsave(filename=f2,plot=read_c,width=10,height=6,dpi=300)
paste0("Generated ", f2)

###### 3. Number of anchors in subregion, colored by read coverage (binned)
anchor_stats_with_ref_paths_only <- anchor_stats_with_ref_paths_only %>%
  mutate(Read_count_bin = factor(
    cut(bp_matched_reads, 
        breaks = c(-Inf, 0, 3, 10, Inf), 
        labels = c("0", "1-3", "4-10", "10+"), 
        right = TRUE),
    levels = c("10+", "4-10", "1-3", "0")  # Set stacking order explicitly
  ))

ar <- ggplot(anchor_stats_with_ref_paths_only, aes(x = Anchor_pos_in_ref_path, fill = Read_count_bin)) +
  geom_histogram(binwidth = 5000, position = "stack", color='black') +  # Dynamic binning
  labs(x = paste("Anchor position in", chr, "reference path"), y = "Number of Anchors", fill = "Read Coverage") +
  theme_bw() +
  scale_x_continuous(breaks =seq(start-5000,end+5000,5000)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) +
  scale_fill_brewer(palette = "Set2")

f3 <- paste0(output_prefix, "_anchor-readcov-dist.png")
ggsave(filename=f3,plot=ar,width=10,height=6,dpi=300)
paste0("Generated ", f3)

cat("anchor to read coverage stats summarised")
anchor_stats_with_ref_paths_only %>% group_by(Read_count_bin) %>% summarise(count=n())

###### 4. Output a BED file with coordinates of anchors to load on UCSC-GB
bed_df <- anchor_stats_with_ref_paths_only[,c('Anchor_path', 'Anchor_length', 'chrom', 'start', 'end', 'Read_count_bin')]

bed_df <- bed_df %>%
    mutate(
        name=paste0(Anchor_path,"_",Anchor_length),
        score=0,
        strand="+",
        thickStart = start,
        thickEnd = end, 
        itemRgb = ifelse(Read_count_bin=="1-3", "141,160,203", ifelse(Read_count_bin=="4-10", "252,141,98", ifelse(Read_count_bin=="0", "231,138,195", "102,194,165"))),
        blockCount = 1, 
        blockSizes = end - start, 
        blockStarts = 0
    ) %>%
    dplyr::select(chrom, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts)

f4 <- paste0(output_prefix, "_anchors_final.bed")
write.table(bed_df, file = f4, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
paste0("Generated ", f4)

####### 5. Add zygosity column to the dataframe and save the master dataframe as TSV
anchor_stats_all <- anchor_stats_all %>% 
    mutate(
    snarl_rank = sapply(
      strsplit(snarl_id, "-"),
      function(parts) {
        nums <- as.numeric(parts[grepl("^\\d+$", parts)])  # keep only numeric strings
        min(nums, na.rm = TRUE)
      }
    )
  ) %>%  # Extract first numeric component
  mutate(snarl_rank=as.integer(snarl_rank)) %>%
  arrange(snarl_rank) %>% group_by(snarl_id) %>% 
  mutate(zygosity = n()) %>% ungroup() %>%
  dplyr::select(Anchor_path, Anchor_length, bp_matched_reads, snarl_id, snarl_rank, everything())
f5 <- paste0(output_prefix, "_anchors_master_table.tsv")
write.table(anchor_stats_all, file = f5, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
paste0("Generated ", f5)

########## Analyze anchor usage by Shasta
# Categorize anchors based on forward_id
anchor_stats_all <- anchor_stats_all %>%
  mutate(
    shasta_status = case_when(
      is.na(forward_id) ~ "NA",
      forward_id == "NA" ~ "NA",
      grepl("Discarded due to coverage", forward_id) ~ "Rejected_low_coverage",
      grepl("Discarded when clipping to markers", forward_id) ~ "Rejected_clipping_to_markers",
      TRUE ~ "Used"
    )
  )

# Count anchors by status
anchor_usage_stats <- anchor_stats_all %>%
  group_by(shasta_status) %>%
  summarise(count = n()) %>%
  ungroup()

# Print usage statistics
cat("\nAnchor Usage Statistics:\n")
print(anchor_usage_stats)

# Save usage statistics to file
f_stats <- paste0(output_prefix, "_anchor_usage_stats.txt")
write.table(anchor_usage_stats, file = f_stats, sep = "\t", row.names = FALSE, quote = FALSE)
paste0("Generated ", f_stats)

# Save a summary TSV with key anchor status columns
anchor_status_summary <- anchor_stats_all %>%
  dplyr::select(Anchor_path, snarl_id, zygosity, Anchor_pos_in_ref_path, shasta_status)
f_summary <- paste0(output_prefix, "_anchor_status_summary.tsv")
write.table(anchor_status_summary, file = f_summary, sep = "\t", row.names = FALSE, quote = FALSE)
paste0("Generated ", f_summary)

# Filter out anchors with no reference path position
anchor_stats_for_plot <- anchor_stats_all %>%
  filter(!is.na(Anchor_pos_in_ref_path) & Anchor_pos_in_ref_path > 0)

# Create a new plot showing anchor positions with usage status
usage_plot <- ggplot(anchor_stats_for_plot, aes(x = Anchor_pos_in_ref_path, y = zygosity, color = shasta_status)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(x = "Anchor position in reference path", 
       y = "Zygosity", 
       title = "Anchor positions colored by Shasta usage status",
       color = "Usage Status") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = sort(unique(anchor_stats_for_plot$zygosity))) +
  scale_x_continuous(breaks = seq(min(anchor_stats_for_plot$Anchor_pos_in_ref_path), 
                                max(anchor_stats_for_plot$Anchor_pos_in_ref_path), 
                                by = 5000)) +
  scale_color_manual(values = c("Used" = "#707070", "Rejected_low_coverage" = "#298c8c", "Rejected_clipping_to_markers" = "#800074", "NA" = "#f1a226"))

f7 <- paste0(output_prefix, "_refPos_zygosity.png")
ggsave(filename = f7, plot = usage_plot, width = 10, height = 6, dpi = 300)
paste0("Generated ", f7)

########## Output some stats
m1.df <- anchor_stats_all %>%
  filter(zygosity == 2) %>%
  arrange(snarl_rank)
m2.df <- m1.df %>% group_by(snarl_id) %>% filter(all(bp_matched_reads >= 3)) %>% ungroup()

outfile <- paste0(output_prefix, "_logs.txt")
cat("Number of anchors given to Shasta: ", nrow(anchor_stats_all), "\n",
    "Number of anchors >= 30 bps: ", sum(anchor_stats_all$Anchor_length >= 30), "\n",
    "Number of het-snarls >= 30 bps: ", sum(m1.df$Anchor_length >= 30) %/% 2, "\n",
    "Number of het-snarls >= 30 bps with >=3 read coverage on both anchors (used in phasing): ", sum(m2.df$Anchor_length >= 30) %/% 2, "\n",
    "Number of anchors >= 100 bps: ", sum(anchor_stats_all$Anchor_length >= 100), "\n",
    "Number of het-snarls >= 100 bps: ", sum(m1.df$Anchor_length >= 100) %/% 2, "\n",
    "Number of het-snarls >= 100 bps with >=3 read coverage on both anchors (used in phasing): ", sum(m2.df$Anchor_length >= 100) %/% 2,
    file = outfile)

####### 6. Snarl zygosity
# Plot snarl_id vs. zygosity, colored by shasta_status
# Use the first anchor's shasta_status for each snarl_id
zygosity_summary_a_col <- anchor_stats_all %>%
  arrange(snarl_rank) %>%
  group_by(snarl_id) %>%
  summarise(zygosity = n(), snarl_rank = first(snarl_rank), shasta_status = first(shasta_status)) %>%
  ungroup() %>%
  arrange(snarl_rank) %>%
  mutate(snarl_id = factor(snarl_id, levels = unique(snarl_id)))

a <- ggplot(zygosity_summary_a_col, aes(x = snarl_rank, y = zygosity, color = shasta_status)) +
  geom_point(alpha=0.6, size=2) +
  labs(x = "Snarl ID", y = "Zygosity", title = paste("Number of anchors per snarl ID", args$region), color = "Usage Status") +
  scale_y_continuous(breaks = seq(1, max(zygosity_summary_a_col$zygosity)+1, 1)) +
  scale_x_continuous(breaks = seq(0, max(zygosity_summary_a_col$snarl_rank)+1, 20)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c("Used" = "#707070", "Rejected_low_coverage" = "#298c8c", "Rejected_clipping_to_markers" = "#800074", "NA" = "#f1a226"))
f6 <- paste0(output_prefix, "_snarlID_zygosity.png")
ggsave(filename=f6,plot=a,width=10,height=6,dpi=300)
paste0("Generated ", f6)

######### 7. Read coverage per snarl (PDF file)
anchors_filtered <- anchor_stats_all %>% 
  filter(!is.na(nodes_matching_reads)) %>%
  arrange(snarl_rank) %>%
  mutate(forward_id = ifelse(is.na(forward_id), paste0("NA"), forward_id))

anchors_filtered$bp_matched_reads[is.na(anchors_filtered$bp_matched_reads)] <- 0

# Split the dataframe into 100 parts
split_data <- split(anchors_filtered, cut(seq_len(nrow(anchors_filtered)), 100, labels = FALSE))

# Open a PDF file to save all plots
f8 <- paste0(output_prefix, "_snarl_reads_stats.pdf")
pdf(f8, width = 12, height = 5)
paste0("Generated ", f8)

# Loop through each part and create plots
for (i in seq_along(split_data)) {
  anchors_long <- split_data[[i]] %>%
    pivot_longer(cols = c(bp_matched_reads), 
                 names_to = "read_type", 
                 values_to = "read_count")

  # Create stacked bar plot with dodged forward_id
  coverage_plot <- ggplot(anchors_long, aes(x = interaction(snarl_id, forward_id), y = read_count, fill = read_type)) +
    geom_bar(stat = "identity", position = "stack") +  # Stack by read_type
    geom_text(aes(label = read_count), 
        position = position_stack(vjust = 0.5), 
        size = 3, color = "black") +
    facet_grid(~snarl_rank, scales = "free_x") +          # Facet by snarl_id for side-by-side bars
    labs(x = "Snarl_ID.Shasta_anchor_ID", y = "Read Counts", 
         title = paste("Reads mapped to each anchor (Part", i, "of 100)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",                       # Move legend to the top
          legend.direction = "horizontal") +             # Arrange legend items horizontally
    scale_fill_manual(values = c("bp_matched_reads" = "#69b3a2"))

  # Print each plot to the PDF
  print(coverage_plot)
}
# Close the PDF file
dev.off()


######### 8. MAPQ distribution plot
read_processed_df <- read.table(read_processed_tsv_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("read_name", "read_len", "relative_strand", "mapq", "div", "path_start", "path_end", "node_list", "orientation_list", "cs_line"))

read_processed_df <- read_processed_df %>% select(read_name, read_len, relative_strand, mapq, div)
message("Number of reads aligned to the region: ", nrow(read_processed_df))

f9 <- paste0(output_prefix, "_mapq_distribution.pdf")

# Plot 1: Proportion of MAPQ=60
p1 <- read_processed_df %>%
  mutate(is_mapq60 = ifelse(mapq == 60, "MAPQ 60", "Other")) %>%
  ggplot(aes(x = is_mapq60)) +
  geom_bar(aes(y = after_stat(count / sum(count))),
           fill = "#428698", alpha = 1) +
  labs(
    x = "",
    y = "Proportion of Reads"
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))

# Plot 2: Distribution excluding MAPQ=60
p2 <- read_processed_df %>%
  filter(mapq != 60) %>%
  ggplot(aes(x = mapq)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, fill = "#428698", alpha = 1
  ) +
  labs(
    x = "",
    y = ""
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))


# Combine side by side
(p1 | p2) + 
  plot_annotation(
    title = "MAPQ Distribution",
    theme = theme(plot.title = element_text(face = "bold", size = 22, hjust = 0.5))
  )

# Save the combined plot
ggsave(filename=f9,plot=p1 + p2,width=8,height=6,dpi=400)
paste0("Generated ", f9)
