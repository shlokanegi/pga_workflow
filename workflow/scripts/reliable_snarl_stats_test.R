#!/usr/bin/env Rscript

# This script takes in the reliable snarls TSV file and the snarl compatibility JSON file
# and generates the following plots and tables:
# 1. Snarl compatibility plot
# 2. Reliable snarl counts
# 3. Reliable snarl zygosity counts
# 4. Snarl compatibility fractions
# 5. Snarl compatibility fractions plot

library(dplyr)
library(data.table)
library(ggplot2)
library(jsonlite)
library(purrr)
library(tidyr)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript reliable_snarl_stats.R <reliable_snarls_tsv> <snarl_compatibility_json> <snarl_coverage_json> <snarl_allelic_coverage_json> <snarl_allelic_coverage_extended_json> <snarl_variant_types_json> <snarl_positions_tsv> <output_dir>", call. = FALSE)
}

reliable_snarls_tsv <- args[1]
snarl_compatibility_json <- args[2]
snarl_coverage_json <- args[3]
snarl_allelic_coverage_json <- args[4]
snarl_allelic_coverage_extended_json <- args[5]
snarl_variant_types_json <- args[6]
snarl_positions_tsv <- args[7]
output_dir <- args[8]


#### --- Output file/plots filenames --- ####
snarl_compatibility_plot_pdf <- file.path(output_dir, "snarl_compatibility_plot.pdf") # top of x-axis is the snarl compatibility bars and bottom is the coverage bars
snarl_compatibility_plot_refpos_pdf <- file.path(output_dir, "snarl_compatibility_plot_refpos.pdf")
snarl_compatibility_false_stratifications_plot_pdf <- file.path(output_dir, "snarl_compatibility_false_stratifications_plot.pdf")
reliable_snarls_counts_tsv <- file.path(output_dir, "reliable_snarls_counts.tsv")
reliable_snarls_zygosity_counts_tsv <- file.path(output_dir, "reliable_snarls_zygosity_counts.tsv")
snarl_compatibility_fractions_tsv <- file.path(output_dir, "snarl_compatibility_fractions.tsv")
snarl_compatibility_fractions_plot_pdf <- file.path(output_dir, "snarl_compatibility_fractions_plot.pdf")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- Load and Process Data ---
message("Loading reliable snarls data from: ", reliable_snarls_tsv)
reliable_snarls.df <- read.table(reliable_snarls_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
reliable_snarls.df <- reliable_snarls.df %>% dplyr::select(snarl_id, zygosity, is_reliable)

message("Loading snarl compatibility data from: ", snarl_compatibility_json)
snarl_compatibility.list <- jsonlite::fromJSON(snarl_compatibility_json)

message("Loading snarl variant types data from: ", snarl_variant_types_json)
snarl_variant_types.list <- jsonlite::fromJSON(snarl_variant_types_json)

# Snarl coverage data before extension
message("Loading snarl coverage data from: ", snarl_coverage_json)
snarl_coverage.list <- jsonlite::fromJSON(snarl_coverage_json)
message("Loading snarl allelic coverage data from: ", snarl_allelic_coverage_json)
snarl_allelic_coverage.list <- jsonlite::fromJSON(snarl_allelic_coverage_json)
# Snarl coverage data after extension
message("Loading snarl allelic coverage data from: ", snarl_allelic_coverage_extended_json)
snarl_allelic_coverage_extended.list <- jsonlite::fromJSON(snarl_allelic_coverage_extended_json)

message("Loading snarl reference positions from: ", snarl_positions_tsv)
snarl_refpos.df <- read.table(snarl_positions_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Check if required columns exist before proceeding
if (!"snarl_id" %in% names(snarl_refpos.df) || !"Anchor_pos_in_ref_path" %in% names(snarl_refpos.df)) {
  stop("The snarl positions TSV file must contain 'snarl_id' and 'Anchor_pos_in_ref_path' columns.")
}

snarl_refpos.df <- snarl_refpos.df %>%
  rename(ref_pos = Anchor_pos_in_ref_path) %>%
  select(snarl_id, ref_pos)


# Tidy the snarl compatibility data into a long-format data frame.
# This single data frame will be used for all subsequent analysis.
compatibility_long.df <- purrr::map_dfr(names(snarl_compatibility.list), ~{
  vals <- unlist(snarl_compatibility.list[[.x]])
  vals <- vals[!is.na(vals)] # Remove NAs
  
  if (length(vals) == 0) {
    return(NULL) # Skip if no valid data
  }
  
  # Ensure compatibility status is a lowercase character vector and clean up
  vals_chr <- tolower(as.character(vals))
  # Consolidate detailed failure reasons into broader categories by removing extra details
  vals_chr <- sub(" .*$", "", vals_chr)
  
  data.frame(
    source_snarl_id = as.integer(.x),
    is_compatible = vals_chr,
    stringsAsFactors = FALSE
  )
})

# Tidy the snarl variant types data into a long-format data frame.
snarl_variant_types.df <- purrr::map_dfr(names(snarl_variant_types.list), ~{
  vals <- snarl_variant_types.list[[.x]] # SNP, INDEL or MNP
  vals_chr <- tolower(as.character(vals)) # snp, indel or mnp

  data.frame(
    snarl_id = as.integer(.x),
    variant_type = vals_chr,
    stringsAsFactors = FALSE
  )
})

head(snarl_variant_types.df)

# Tidy the snarl coverage data into a long-format data frame.
snarl_coverage.df <- purrr::map_dfr(
  .x = names(snarl_coverage.list),
  .f = ~ tibble::tibble(
    snarl_id = as.integer(.x),
    coverage = snarl_coverage.list[[.x]]
  )
)

allelic_coverage.df <- purrr::map_dfr(names(snarl_allelic_coverage.list), ~{
    snarl_data <- snarl_allelic_coverage.list[[.x]]
    data.frame(
        snarl_id = as.integer(.x),
        allele_id = names(snarl_data),
        coverage = as.integer(unlist(snarl_data)),
        stringsAsFactors = FALSE
    )
})

# Tidy the extended allelic coverage data, handling merged and missing snarls
allelic_coverage_extended.df <- purrr::map_dfr(names(snarl_allelic_coverage_extended.list), ~{
    current_snarl_id_str <- .x
    snarl_data <- snarl_allelic_coverage_extended.list[[current_snarl_id_str]]

    # Split merged snarl IDs (e.g., "1-2") into individual integer IDs
    original_snarl_ids <- as.integer(unlist(strsplit(current_snarl_id_str, "-")))

    # Create a dataframe for the allele data for the current snarl
    alleles_df <- data.frame(
        allele_id = names(snarl_data),
        coverage = as.integer(unlist(snarl_data)),
        stringsAsFactors = FALSE
    )

    # Create rows for each original snarl ID, duplicating allele/coverage info
    tidyr::crossing(snarl_id = original_snarl_ids, alleles_df)
})

# Identify snarls present before extension but missing after
all_snarl_ids_before_extension <- unique(snarl_coverage.df$snarl_id)
all_snarl_ids_after_extension <- unique(allelic_coverage_extended.df$snarl_id)
missing_snarl_ids <- setdiff(all_snarl_ids_before_extension, all_snarl_ids_after_extension)

# Add missing snarls with coverage=0 for allele 0
if (length(missing_snarl_ids) > 0) {
    missing_snarls.df <- data.frame(
        snarl_id = missing_snarl_ids,
        allele_id = "0",
        coverage = 0
    )
    allelic_coverage_extended.df <- dplyr::bind_rows(allelic_coverage_extended.df, missing_snarls.df)
}

# Get all unique allele IDs from both dataframes to create a shared legend
all_allele_ids <- unique(c(
  as.character(allelic_coverage.df$allele_id),
  as.character(allelic_coverage_extended.df$allele_id)
))
# Sort numerically to ensure consistent ordering
all_allele_ids <- all_allele_ids[order(as.numeric(all_allele_ids))]

# Apply the same factor levels to both dataframes to ensure a single, combined legend
allelic_coverage.df$allele_id <- factor(allelic_coverage.df$allele_id, levels = all_allele_ids)
allelic_coverage_extended.df$allele_id <- factor(allelic_coverage_extended.df$allele_id, levels = all_allele_ids)

#########################################################################################################
################################### PLOTS and STATS #####################################################
#########################################################################################################

#### --- 1. Snarl compatibility analysis --- ####
# Create a summary of compatibility counts for the first plot
compatibility_counts.df <- compatibility_long.df %>%
  group_by(source_snarl_id, is_compatible) %>%
  summarise(count = n(), .groups = 'drop') %>%
  rename(snarl_id = source_snarl_id)

# Calculate total counts for positioning the points on top of the bars
total_counts.df <- compatibility_counts.df %>%
  group_by(snarl_id) %>%
  summarize(total_count = sum(count), .groups = 'drop')

# Join with variant types
variant_points.df <- snarl_variant_types.df %>%
  left_join(total_counts.df, by = "snarl_id")

# Calculate the total coverage per snarl for both allelic coverage dataframes
total_coverage_pre <- allelic_coverage.df %>%
  group_by(snarl_id) %>%
  summarise(total_coverage = sum(coverage))

total_coverage_post <- allelic_coverage_extended.df %>%
  group_by(snarl_id) %>%
  summarise(total_coverage = sum(coverage))

# Find the overall maximum coverage to synchronize the y-axes of the coverage plots
max_coverage_snarl_id <- max(max(total_coverage_pre$total_coverage), max(total_coverage_post$total_coverage), na.rm = TRUE)


# --- 1st plot: Snarl Compatibility (Aggregated) ---
# Aggregate all 'false' categories into one for the main plot
compatibility_counts_aggregated.df <- compatibility_long.df %>%
  mutate(is_compatible = ifelse(startsWith(is_compatible, "false"), "false", is_compatible)) %>%
  group_by(source_snarl_id, is_compatible) %>%
  summarise(count = n(), .groups = 'drop') %>%
  rename(snarl_id = source_snarl_id)

# Recalculate total counts for positioning points on the aggregated plot
total_counts_aggregated.df <- compatibility_counts_aggregated.df %>%
  group_by(snarl_id) %>%
  summarize(total_count = sum(count), .groups = 'drop')

# Join with variant types for the aggregated plot
variant_points_aggregated.df <- snarl_variant_types.df %>%
  left_join(total_counts_aggregated.df, by = "snarl_id")

p_compatibility <- compatibility_counts_aggregated.df %>%
  ggplot(aes(x = snarl_id, y = count, fill = is_compatible)) +
  geom_col(position = "stack", width = 1) +
  geom_point(data = variant_points_aggregated.df, 
             aes(x = snarl_id, y = total_count, shape = variant_type, fill = NULL), 
             color = "black", size = 2.0) +
  geom_point(data = variant_points_aggregated.df, 
             aes(x = snarl_id, y = total_count, color = variant_type, shape = variant_type, fill = NULL), 
             size = 1.5, alpha = 1) +
  scale_fill_manual(
    name = "Snarl Compatibility Status",
    values = c("true" = "#82bdd0", "false" = "#e15759"),
    labels = c("true" = "Compatible", "false" = "Incompatible")
  ) +
  labs(
    y = "Number of Linked Snarls",
    x = NULL,  # Remove x-axis title for this plot
    color = "Variant Type",
    shape = "Variant Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 5.5, 1, 5.5)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(compatibility_counts_aggregated.df$snarl_id), by = ceiling(max(compatibility_counts_aggregated.df$snarl_id)/15)),
    labels = scales::comma
  ) +
  scale_y_continuous(
    breaks = seq(0, max(variant_points_aggregated.df$total_count, na.rm = TRUE) + 20, by = ceiling(max(variant_points_aggregated.df$total_count, na.rm = TRUE) / 10)),
    labels = scales::comma,
    expand = c(0, 1.5)
  ) +
  scale_color_manual(values = c("snp" = "#bc272d", "mnp" = "#50a89f", "indel" = "#1630c1")) +
  scale_shape_manual(values = c("snp" = 16, "mnp" = 15, "indel" = 17))

# --- New plot: False_* Stratifications ---
# Filter for false categories only
compatibility_counts_false.df <- compatibility_counts.df %>%
  filter(startsWith(is_compatible, "false"))

# Total counts for positioning points on the false plot
total_counts_false.df <- compatibility_counts_false.df %>%
  group_by(snarl_id) %>%
  summarize(total_count = sum(count), .groups = 'drop')

# Join with variant types for the false plot
variant_points_false.df <- snarl_variant_types.df %>%
  left_join(total_counts_false.df, by = "snarl_id")

# Dynamically create color palette for the false plot
library(RColorBrewer)
all_statuses <- unique(compatibility_counts_false.df$is_compatible)
false_statuses <- sort(all_statuses[startsWith(all_statuses, "false")])

if (length(false_statuses) > 0) {
    num_colors <- length(false_statuses)
    # Get a palette of distinct colors
    if (num_colors <= 9) {
        distinct_colors <- brewer.pal(max(3, num_colors), "Set1")[1:num_colors]
    } else {
        distinct_colors <- colorRampPalette(brewer.pal(9, "Set1"))(num_colors)
    }
    names(distinct_colors) <- false_statuses
} else {
    distinct_colors <- c()
}

status_labels <- c(
  "false" = "Incompatible",
  "false_lowcov" = "Incompatible (lowCov)",
  "false_setsunequal" = "Incompatible (setsUnequal)",
  "false_gtestfailed" = "Incompatible (g-testFailed)",
  "false_gtestproducednohypotheses" = "Incompatible (noHypotheses)",
  "false_besthypothesisnotabijective" = "Incompatible (notBijection)",
  "false_besthypothesisgtoohigh" = "Incompatible (gTooHigh)",
  "false_hypothesesnotwellseparated" = "Incompatible (notSeparated)"
)

p_compatibility_false <- compatibility_counts_false.df %>%
  ggplot(aes(x = snarl_id, y = count, fill = is_compatible)) +
  geom_col(position = "stack", width = 1) +
  geom_point(data = variant_points_false.df, 
             aes(x = snarl_id, y = total_count, shape = variant_type, fill = NULL), 
             color = "black", size = 2.0) +
  geom_point(data = variant_points_false.df, 
             aes(x = snarl_id, y = total_count, color = variant_type, shape = variant_type, fill = NULL), 
             size = 1.5, alpha = 1) +
  scale_fill_manual(
    name = "Snarl Incompatibility Status",
    values = distinct_colors,
    labels = status_labels[names(distinct_colors)]
  ) +
  labs(
    title = "Snarl Incompatibility Breakdown",
    y = "Number of Linked Snarls",
    x = "Source Snarl ID",
    color = "Variant Type",
    shape = "Variant Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_x_continuous(
    breaks = seq(0, max(compatibility_counts_false.df$snarl_id), by = ceiling(max(compatibility_counts_false.df$snarl_id)/15)),
    labels = scales::comma
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = c(0, 1.5)
  ) +
  scale_color_manual(values = c("snp" = "#bc272d", "mnp" = "#50a89f", "indel" = "#1630c1")) +
  scale_shape_manual(values = c("snp" = 16, "mnp" = 15, "indel" = 17))

ggsave(snarl_compatibility_false_stratifications_plot_pdf, plot = p_compatibility_false, width = 12, height = 7, dpi = 400)


#### --- 2nd plot: Snarl Allelic Coverage --- ####
p_allelic_coverage <- allelic_coverage.df %>%
  ggplot(aes(x = snarl_id, y = coverage, fill = allele_id)) +
  geom_col(position = "stack", width = 1) +
  labs(
    y = "Allelic Coverage\n(pre-extension)",
    x = "Source Snarl ID",
    fill = "Allele ID"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(1, 5.5, 1, 5.5)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(allelic_coverage.df$snarl_id), by = ceiling(max(allelic_coverage.df$snarl_id)/15)),
    labels = scales::comma
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, max_coverage_snarl_id * 1.05)) +
  scale_fill_brewer(palette = "Paired")

#### --- 3rd plot: Snarl Allelic Coverage Extended anchors --- ####
p_allelic_coverage_extended <- allelic_coverage_extended.df %>%
  ggplot(aes(x = snarl_id, y = coverage, fill = allele_id)) +
  geom_col(position = "stack", width = 1) +
  labs(
    y = "Allelic Coverage\n(post-extension)",
    x = "Source Snarl ID",
    fill = "Allele ID"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    plot.margin = margin(1, 5.5, 5.5, 5.5)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(allelic_coverage_extended.df$snarl_id), by = ceiling(max(allelic_coverage_extended.df$snarl_id)/15)),
    labels = scales::comma
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, max_coverage_snarl_id * 1.05)) +
  scale_fill_brewer(palette = "Paired")

# --- Combine and Save All Plots using cowplot ---
# Create a title grob for the entire figure
title_grob <- ggdraw() + 
    draw_label("Snarl compatibility and allelic coverage", fontface = 'bold', size = 18, x = 0.5)

# Extract the legends from the plots
leg_compat <- get_legend(
  p_compatibility +
    guides(
      fill = guide_legend(nrow = 2, title.position = "left"),
      color = guide_legend(nrow = 2, title.position = "left"),
      shape = guide_legend(nrow = 2, title.position = "left")
    ) +
    theme(legend.box = "horizontal")
)
leg_allelic <- get_legend(p_allelic_coverage + guides(fill = guide_legend(nrow = 2, title.position = "left")))

# Arrange the legends horizontally
legs <- plot_grid(leg_compat, leg_allelic, ncol = 2, align = 'h', rel_widths = c(2, 1))

# Arrange the four plots vertically, with their legends removed
plots <- plot_grid(
  p_compatibility + theme(legend.position = "none"),
  p_allelic_coverage + theme(legend.position = "none"),
  p_allelic_coverage_extended + theme(legend.position = "none"),
  align = 'v',
  ncol = 1,
  rel_heights = c(1.5, 0.8, 1)
)

# Combine the title and the plots
main_panel <- plot_grid(title_grob, plots, ncol = 1, rel_heights = c(0.08, 1))

# Combine the main plot panel with the legends at the bottom
final_plot <- plot_grid(main_panel, legs, ncol = 1, rel_heights = c(1, .1))

ggsave(snarl_compatibility_plot_pdf, plot = final_plot, width = 12, height = 9, dpi = 400)


# --- Create the new plot with reference positions ---

# --- Data Preparation for all refpos plots ---
# Master list of all snarls with their unique reference positions and zygosity.
master_snarls_with_refpos <- snarl_refpos.df %>%
  distinct(snarl_id, .keep_all = TRUE) %>%
  left_join(reliable_snarls.df, by = "snarl_id")

# Master list of all snarls with reference positions and allelic coverage (pre-extension)
master_snarls_with_refpos_coverage <- master_snarls_with_refpos %>%
  left_join(allelic_coverage.df, by = "snarl_id") %>%
  mutate(coverage = ifelse(is.na(coverage), 0, coverage)) %>%
  filter(!is.na(ref_pos) & ref_pos > 0)

# Master list of all snarls with reference positions and allelic coverage (post-extension)
master_snarls_with_refpos_coverage_extended <- master_snarls_with_refpos %>%
  left_join(allelic_coverage_extended.df, by = "snarl_id") %>%
  mutate(coverage = ifelse(is.na(coverage), 0, coverage)) %>%
  filter(!is.na(ref_pos) & ref_pos > 0)

# --- Recalculate Max Coverage for Y-axis scaling ---
total_coverage_pre_full <- master_snarls_with_refpos_coverage %>%
  group_by(snarl_id) %>%
  summarise(total_coverage = sum(coverage))

total_coverage_post_full <- master_snarls_with_refpos_coverage_extended %>%
  group_by(snarl_id) %>%
  summarise(total_coverage = sum(coverage))

max_coverage_refpos <- max(max(total_coverage_pre_full$total_coverage, na.rm = TRUE), max(total_coverage_post_full$total_coverage, na.rm = TRUE), na.rm = TRUE)

# --- Plotting Section ---
# The data for the top plot (compatibility) is prepared here. It correctly uses only het snarls.
snarl_refpos_unique <- snarl_refpos.df %>%
  distinct(snarl_id, .keep_all = TRUE)
  
refpos_data <- compatibility_long.df %>%
  left_join(snarl_refpos_unique, by = c("source_snarl_id" = "snarl_id")) %>%
  rename(snarl_id = source_snarl_id)
  
compatibility_counts_aggregated_refpos.df <- refpos_data %>%
    mutate(is_compatible = ifelse(startsWith(is_compatible, "false"), "false", is_compatible)) %>%
    filter(!is.na(ref_pos) & ref_pos > 0) %>%
    group_by(snarl_id, is_compatible, ref_pos) %>%
    summarise(count = n(), .groups = 'drop')

total_counts_aggregated_refpos.df <- compatibility_counts_aggregated_refpos.df %>%
    group_by(snarl_id, ref_pos) %>%
    summarize(total_count = sum(count), .groups = 'drop')

variant_points_aggregated_refpos.df <- snarl_variant_types.df %>%
    right_join(total_counts_aggregated_refpos.df, by = "snarl_id")

plot_data_refpos <- compatibility_counts_aggregated_refpos.df %>%
  filter(!is.na(ref_pos) & ref_pos > 0)

# Calculate x-axis limits and tick span
min_pos <- min(plot_data_refpos$ref_pos, na.rm = TRUE)
max_pos <- max(plot_data_refpos$ref_pos, na.rm = TRUE)
tick_span <- round((max_pos - min_pos) / 15)

# --- 1st plot: Snarl Compatibility (refpos, het only) ---
p_compatibility_refpos <- plot_data_refpos %>%
  ggplot(aes(x = ref_pos, y = count, fill = is_compatible)) +
  geom_col(position = "stack", width = 500) +
  geom_point(data = variant_points_aggregated_refpos.df %>% filter(!is.na(total_count)),
             aes(x = ref_pos, y = total_count, shape = variant_type, fill = NULL),
             color = "black", size = 2.0) +
  geom_point(data = variant_points_aggregated_refpos.df %>% filter(!is.na(total_count)),
             aes(x = ref_pos, y = total_count, color = variant_type, shape = variant_type, fill = NULL),
             size = 1.5, alpha = 1) +
  scale_fill_manual(
    name = "Snarl Compatibility Status",
    values = c("true" = "#82bdd0", "false" = "#e15759"),
    labels = c("true" = "Compatible", "false" = "Incompatible")
  ) +
  labs(
    y = "Number of Linked Snarls",
    x = NULL, # No x-axis title for this plot
    color = "Variant Type",
    shape = "Variant Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(), # No x-axis labels
    axis.ticks.x = element_blank(), # No x-axis ticks
    plot.margin = margin(0, 5.5, 1, 5.5)
  ) +
  scale_x_continuous(labels = function(x) scales::comma(x/1000, accuracy = 0.01), breaks = seq(from = min_pos, to = max_pos, by = tick_span), expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = c("snp" = "#bc272d", "mnp" = "#50a89f", "indel" = "#1630c1")) +
  scale_shape_manual(values = c("snp" = 16, "mnp" = 15, "indel" = 17))

# --- 2nd plot: Snarl Allelic Coverage (refpos, all snarls) ---
p_allelic_coverage_refpos <- master_snarls_with_refpos_coverage %>%
  ggplot(aes(x = ref_pos, y = coverage, fill = allele_id)) +
  geom_col(position = "stack", width = 500) +
  labs(
    y = "Allelic Coverage\n(pre-extension)",
    x = NULL, # No x-axis title
    fill = "Allele ID"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(), # No x-axis labels
    axis.ticks.x = element_blank(), # No x-axis ticks
    plot.margin = margin(1, 5.5, 1, 5.5)
  ) +
  scale_x_continuous(labels = function(x) scales::comma(x/1000, accuracy = 0.01), breaks = seq(from = min_pos, to = max_pos, by = tick_span), expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, max_coverage_refpos * 1.05)) +
  scale_fill_brewer(palette = "Paired")

# --- 3rd plot: Snarl Allelic Coverage Extended anchors (refpos, all snarls) ---
p_allelic_coverage_extended_refpos <- master_snarls_with_refpos_coverage_extended %>%
  ggplot(aes(x = ref_pos, y = coverage, fill = allele_id)) +
  geom_col(position = "stack", width = 500) +
  labs(
    y = "Allelic Coverage\n(post-extension)",
    x = "Reference Position (CHM13) in Kbp",
    fill = "Allele ID"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust=1),
    plot.margin = margin(1, 5.5, 5.5, 5.5)
  ) +
  scale_x_continuous(labels = function(x) scales::comma(x/1000, accuracy = 0.01), breaks = seq(from = min_pos, to = max_pos, by = tick_span), expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, max_coverage_refpos * 1.05)) +
  scale_fill_brewer(palette = "Paired")

# --- Combine and Save All Plots using cowplot ---
# Create a title grob for the entire figure
title_grob_refpos <- ggdraw() + 
    draw_label("Snarl compatibility and allelic coverage", fontface = 'bold', size = 18, x = 0.5) +
    draw_label(paste("Region Size:", sprintf("%.1f", (max_pos - min_pos) / 1000), " Kbp  |  x-axis tick interval:", sprintf("%.1f", tick_span / 1000), " Kbp"), size = 13, x = 0.5, y=0.2)

# Extract the legends from the plots
leg_compat_refpos <- get_legend(
  p_compatibility_refpos +
    guides(
      fill = guide_legend(nrow = 2, title.position = "left"),
      color = guide_legend(nrow = 2, title.position = "left"),
      shape = guide_legend(nrow = 2, title.position = "left")
    ) +
    theme(legend.box = "horizontal")
)
leg_allelic_refpos <- get_legend(p_allelic_coverage_refpos + guides(fill = guide_legend(nrow = 2, title.position = "left")))

# Arrange the legends horizontally
legs_refpos <- plot_grid(leg_compat_refpos, leg_allelic_refpos, ncol = 2, align = 'h', rel_widths = c(2, 1))

# Arrange the four plots vertically, with their legends removed
plots_refpos <- plot_grid(
  p_compatibility_refpos + theme(legend.position = "none"),
  p_allelic_coverage_refpos + theme(legend.position = "none"),
  p_allelic_coverage_extended_refpos + theme(legend.position = "none"),
  align = 'v',
  ncol = 1,
  rel_heights = c(1, 0.5, 0.8)
)

# Combine the title and the plots
main_panel_refpos <- plot_grid(title_grob_refpos, plots_refpos, ncol = 1, rel_heights = c(0.1, 1))

# Combine the main plot panel with the legends at the bottom
final_plot_refpos <- plot_grid(main_panel_refpos, legs_refpos, ncol = 1, rel_heights = c(1, .1))

ggsave(snarl_compatibility_plot_refpos_pdf, plot = final_plot_refpos, width = 12, height = 9, dpi = 400)


#### --- 3. Reliable snarl TSV --- ####
reliable_snarls.df %>% filter(zygosity!=1) %>% group_by(is_reliable) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(percentage = count / sum(count) * 100) %>%
  write.table(file = reliable_snarls_counts_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

reliable_snarls.df %>% filter(zygosity!=1) %>% group_by(is_reliable, zygosity) %>%
  summarise(count = n(), .groups = 'drop') %>%
  write.table(file = reliable_snarls_zygosity_counts_tsv, sep = "\t", row.names = FALSE, quote = FALSE)


#### --- 4. Snarl compatibility fractions --- ####
# Compute counts and fraction: True / (True + False).
# This uses the `compatibility_long.df` created earlier.
compatibility_fractions.df <- compatibility_long.df %>%
  group_by(source_snarl_id, is_compatible) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(
    names_from = is_compatible,
    values_from = count,
    values_fill = 0
  ) %>%
  rename_with(~ paste0("n_", .), -source_snarl_id)

# Ensure the 'n_true' column exists, even if no 'true' values were found
if (!"n_true" %in% names(compatibility_fractions.df)) {
  compatibility_fractions.df$n_true <- 0
}

# Calculate total counts for fraction calculation
total_counts_for_fraction <- rowSums(select(compatibility_fractions.df, -source_snarl_id), na.rm = TRUE)

# Calculate fraction and add it to the data frame
compatibility_fractions.df <- compatibility_fractions.df %>%
  mutate(
    fraction_compatible = ifelse(total_counts_for_fraction > 0,
                                 n_true / total_counts_for_fraction,
                                 NA_real_)
  )


# Merge with reliable snarl data frame to carry over the zygosity column
compatibility_fractions.df <- compatibility_fractions.df %>%
  left_join(reliable_snarls.df, by = c("source_snarl_id" = "snarl_id")) %>%
  mutate(zygosity = ifelse(is.na(zygosity), 0, zygosity))  # Fill NA zygosity with 0

# Reorder columns to have zygosity and fraction_compatible after main counts
cols_order <- c(
    "source_snarl_id", "zygosity", "n_true",
    sort(names(compatibility_fractions.df)[startsWith(names(compatibility_fractions.df), "n_false")]),
    "fraction_compatible"
)
# Ensure all expected columns are present before reordering
cols_order <- cols_order[cols_order %in% names(compatibility_fractions.df)]
compatibility_fractions.df <- compatibility_fractions.df %>%
  select(all_of(cols_order))


write.table(compatibility_fractions.df,
            file = snarl_compatibility_fractions_tsv,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


#### --- 5. Snarl compatibility fractions plot --- ####
compatibility_fractions_plot.df <- compatibility_fractions.df %>%
  filter(!is.na(fraction_compatible))

p2 <- ggplot(compatibility_fractions_plot.df, aes(x = fraction_compatible)) +
  geom_histogram(
    bins = 40,
    fill = "#4e79a7",
    color = "white",
    alpha = 0.8
  ) +
  labs(
    title = "Distribution of Snarl Compatibility Fractions",
    subtitle = paste0(
      "Each data point = fraction of compatible links for a snarl."
    ),
    x = "Fraction of Compatible Linked Snarls",
    y = "Number of Snarls"
  ) +
  theme_bw(base_size = 14) +
  scale_y_continuous(
    labels = scales::comma
  )
ggsave(snarl_compatibility_fractions_plot_pdf, plot = p2, width = 10, height = 6, dpi = 400)
