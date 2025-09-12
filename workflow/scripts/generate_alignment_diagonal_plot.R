#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(scales)
  library(gridExtra)
  library(ggrepel)
  library(grDevices)
})

# Function to parse PAF file
parse_paf <- function(paf_file) {
  cat("Parsing PAF file:", paf_file, "\n")
  
  # Read PAF file line by line to handle variable column counts
  lines <- readLines(paf_file)
  cat("Read", length(lines), "lines from PAF file\n")
  
  # Parse each line
  max_cols <- 12  # Start with minimum required columns
  
  for (i in seq_along(lines)) {
    if (lines[i] != "" && !startsWith(lines[i], "#")) {
      fields <- strsplit(lines[i], "\t")[[1]]
      max_cols <- max(max_cols, length(fields))
    }
  }
  
  cat("Maximum columns found:", max_cols, "\n")
  
  # Create empty matrix with max columns
  paf_matrix <- matrix(NA, nrow = length(lines), ncol = max_cols, dimnames = list(NULL, paste0("V", 1:max_cols)))
  
  valid_rows <- 0
  for (i in seq_along(lines)) {
    if (lines[i] != "" && !startsWith(lines[i], "#")) {
      fields <- strsplit(lines[i], "\t")[[1]]
      if (length(fields) >= 12) {  # Must have at least 12 standard PAF columns
        valid_rows <- valid_rows + 1
        paf_matrix[valid_rows, 1:length(fields)] <- fields
      }
    }
  }
  
  # Remove empty rows
  paf_matrix <- paf_matrix[1:valid_rows, , drop = FALSE]
  
  # Convert to data frame
  paf_data <- as.data.frame(paf_matrix, stringsAsFactors = FALSE)
  
  # Standard PAF columns
  colnames(paf_data)[1:12] <- c(
    "query_name", "query_length", "query_start", "query_end", "strand",
    "target_name", "target_length", "target_start", "target_end",
    "matches", "alignment_length", "mapping_quality"
  )
  
  # Convert numeric columns
  numeric_cols <- c("query_length", "query_start", "query_end", 
                   "target_length", "target_start", "target_end",
                   "matches", "alignment_length", "mapping_quality")
  
  for (col in numeric_cols) {
    paf_data[[col]] <- as.numeric(paf_data[[col]])
  }
  
  # Name additional tag columns
  if (ncol(paf_data) > 12) {
    for (i in 13:ncol(paf_data)) {
      colnames(paf_data)[i] <- paste0("tag_", i-12)
    }
  }
  
  # Filter for primary alignments only (check for tp:A:P tag or absence of secondary/supplementary flags)
  primary_indices <- rep(TRUE, nrow(paf_data))
  
  # Check for tp:A: tags (P=primary, S=secondary, I=inversion)
  for (i in 1:nrow(paf_data)) {
    row_tags <- as.character(paf_data[i, 13:ncol(paf_data)])
    tp_tags <- row_tags[grepl("^tp:A:", row_tags)]
    if (length(tp_tags) > 0) {
      tp_value <- gsub("^tp:A:", "", tp_tags[1])
      if (tp_value != "P") {
        primary_indices[i] <- FALSE
      }
    }
  }
  
  paf_data <- paf_data[primary_indices, ]
  
  cat("Found", nrow(paf_data), "valid primary alignments\n")
  return(paf_data)
}

# Function to create detailed alignment segments from CS string
create_alignment_segments <- function(cs_string, target_start, query_start, target_end, query_end) {
  segments <- data.frame()
  variants <- data.frame()
  
  if (is.na(cs_string) || cs_string == "") {
    # No CS string - create simple segment
    segments <- data.frame(
      target_start = target_start,
      target_end = target_end,
      query_start = query_start,
      query_end = query_end
    )
    return(list(segments = segments, variants = variants))
  }
  
  # Remove cs:Z: prefix if present
  cs_string <- gsub("^cs:Z:", "", cs_string)
  
  target_pos <- target_start
  query_pos <- query_start
  segment_start_target <- target_start
  segment_start_query <- query_start
  
  # Parse CS string operations
  operations <- str_extract_all(cs_string, "(:[0-9]+|\\*[acgtACGT][acgtACGT]|\\+[acgtACGT]+|\\-[acgtACGT]+)")[[1]]
  
  for (op in operations) {
    if (str_detect(op, "^:")) {
      # Match - continue segment
      match_len <- as.numeric(str_extract(op, "[0-9]+"))
      target_pos <- target_pos + match_len
      query_pos <- query_pos + match_len
      
    } else if (str_detect(op, "^\\*")) {
      # Substitution/mismatch - end current segment and start new one
      if (target_pos > segment_start_target) {
        segments <- rbind(segments, data.frame(
          target_start = segment_start_target,
          target_end = target_pos,
          query_start = segment_start_query,
          query_end = query_pos
        ))
      }
      
      # Add variant
      variants <- rbind(variants, data.frame(
        target_pos = target_pos,
        query_pos = query_pos,
        type = "mismatch",
        ref_base = str_sub(op, 2, 2),
        alt_base = str_sub(op, 3, 3)
      ))
      
      target_pos <- target_pos + 1
      query_pos <- query_pos + 1
      segment_start_target <- target_pos
      segment_start_query <- query_pos
      
    } else if (str_detect(op, "^\\+")) {
      # Insertion in query - end current segment
      if (target_pos > segment_start_target) {
        segments <- rbind(segments, data.frame(
          target_start = segment_start_target,
          target_end = target_pos,
          query_start = segment_start_query,
          query_end = query_pos
        ))
      }
      
      ins_seq <- str_extract(op, "[acgtACGT]+")
      ins_len <- nchar(ins_seq)
      
      # Add variant
      variants <- rbind(variants, data.frame(
        target_pos = target_pos,
        query_pos = query_pos,
        type = "insertion",
        ref_base = "",
        alt_base = ins_seq
      ))
      
      query_pos <- query_pos + ins_len
      segment_start_target <- target_pos
      segment_start_query <- query_pos
      
    } else if (str_detect(op, "^\\-")) {
      # Deletion in query - end current segment
      if (target_pos > segment_start_target) {
        segments <- rbind(segments, data.frame(
          target_start = segment_start_target,
          target_end = target_pos,
          query_start = segment_start_query,
          query_end = query_pos
        ))
      }
      
      del_seq <- str_extract(op, "[acgtACGT]+")
      del_len <- nchar(del_seq)
      
      # Add variant
      variants <- rbind(variants, data.frame(
        target_pos = target_pos,
        query_pos = query_pos,
        type = "deletion",
        ref_base = del_seq,
        alt_base = ""
      ))
      
      target_pos <- target_pos + del_len
      segment_start_target <- target_pos
      segment_start_query <- query_pos
    }
  }
  
  # Add final segment
  if (target_pos > segment_start_target) {
    segments <- rbind(segments, data.frame(
      target_start = segment_start_target,
      target_end = target_pos,
      query_start = segment_start_query,
      query_end = query_pos
    ))
  }
  
  return(list(segments = segments, variants = variants))
}

# Function to create diagonal plot for a single target contig (returns plot object only)
create_diagonal_plot <- function(alignments, target_contig) {
  cat("Creating plot for target contig:", target_contig, "\n")
  
  # Filter alignments for this target contig
  target_alignments <- alignments[alignments$target_name == target_contig, ]
  
  if (nrow(target_alignments) == 0) {
    cat("No alignments found for", target_contig, "\n")
    return(NULL)
  }
  
  # Create detailed plot data with proper alignment segments
  plot_data <- data.frame()
  variant_data <- data.frame()
  label_data <- data.frame()
  arrow_data <- data.frame()  # separate data for arrows
  box_data <- data.frame()    # NEW: data for grey boxes around alignments
  
  for (i in 1:nrow(target_alignments)) {
    aln <- target_alignments[i, ]
    
    # Find CS string for detailed alignment
    cs_col <- which(sapply(aln, function(x) !is.na(x) && grepl("^cs:Z:", x)))
    cs_string <- ifelse(length(cs_col) > 0, as.character(aln[cs_col[1]]), NA)
    
    # Create detailed segments
    result <- create_alignment_segments(cs_string, aln$target_start, aln$query_start, aln$target_end, aln$query_end)
    
    # Add segments (without arrows)
    if (nrow(result$segments) > 0) {
      segments <- result$segments
      segments$query_name <- aln$query_name
      segments$strand <- aln$strand
      plot_data <- rbind(plot_data, segments)
    }
    
    # Add variants
    if (nrow(result$variants) > 0) {
      variants <- result$variants
      variants$query_name <- aln$query_name
      variant_data <- rbind(variant_data, variants)
    }
    
    # Add label data (one label per alignment, positioned at the middle)
    mid_target <- (aln$target_start + aln$target_end) / 2
    mid_query <- (aln$query_start + aln$query_end) / 2
    label_data <- rbind(label_data, data.frame(
      target_pos = mid_target,
      query_pos = mid_query,
      query_name = aln$query_name
    ))
    
    # Add arrow data (one arrow per alignment, at the end)
    arrow_data <- rbind(arrow_data, data.frame(
      query_start = aln$query_end - (aln$query_end - aln$query_start) * 0.02,  # Start arrow 2% before end
      target_start = aln$target_end - (aln$target_end - aln$target_start) * 0.02,
      query_end = aln$query_end,
      target_end = aln$target_end,
      query_name = aln$query_name
    ))
    
    # Add box data for each alignment region (exact fit to alignment)
    box_data <- rbind(box_data, data.frame(
      xmin = aln$query_start,
      xmax = aln$query_end,
      ymin = aln$target_start,
      ymax = aln$target_end,
      query_name = aln$query_name
    ))
  }
  
  # Get target contig length
  target_length <- max(target_alignments$target_length)
  
  # Create the plot
  p <- ggplot() +
    # Add light grey boxes first (so they appear behind everything else)
    geom_rect(data = box_data,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "lightgrey", alpha = 0.3, color = "lightgrey", linewidth = 0.2) +
    # Add alignment segments without arrows (all in same color)
    geom_segment(data = plot_data, 
                 aes(x = query_start, y = target_start, 
                     xend = query_end, yend = target_end),
                 color = "#7d8091", linewidth = 1, alpha = 0.7) +
    # Add arrows at the end of each alignment
    geom_segment(data = arrow_data,
                 aes(x = query_start, y = target_start,
                     xend = query_end, yend = target_end),
                 color = "#7d8091", linewidth = 1, alpha = 0.7,
                 arrow = arrow(length = unit(0.15, "inches"), type = "closed"))
  
  # Add variant points and scales only if variants exist
  if (nrow(variant_data) > 0) {
    p <- p +
      geom_point(data = variant_data,
                 aes(x = query_pos, y = target_pos, color = type),
                 size = 2.5, alpha = 0.8, stroke = 0) +
      scale_color_manual(values = c("mismatch" = "#E41A1C", "insertion" = "#377EB8", "deletion" = "#4DAF4A"),
                         name = "Variant Type")
  }
  
  # Add contig name labels with repel (in boxes, away from diagonal)
  if (nrow(label_data) > 0) {
    p <- p +
      geom_label_repel(data = label_data,
                       aes(x = query_pos, y = target_pos, label = query_name),
                       size = 3, fontface = "bold", 
                       box.padding = 0.5, point.padding = 0.3,
                       max.overlaps = Inf, force = 3,
                       nudge_y = target_length * 0.05,  # Nudge labels away from diagonal
                       fill = "white", alpha = 0.8,
                       label.padding = unit(0.2, "lines"),
                       label.r = unit(0.15, "lines"))
  }
  
  # Add remaining formatting
  p <- p +
    labs(
      title = paste("Alignment Plot:", target_contig),
      subtitle = paste("Target length:", format(target_length, big.mark = ","), "| Primary alignments only | Grey boxes show alignment regions"),
      x = "Position in Shasta Assembly Contig (bp)",
      y = "Position in Hifiasm Assembly Contig (bp)"
    ) +
    
    # Format axes with commas
    scale_x_continuous(labels = comma_format()) +
    scale_y_continuous(labels = comma_format()) +
    
    # Theme
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Main function
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 3) {
    cat("Usage: Rscript generate_alignment_diagonal_plot.R <paf_file> <output_dir> <sample_prefix>\n")
    cat("Arguments:\n")
    cat("  paf_file: Path to PAF alignment file\n")
    cat("  output_dir: Output directory for plots\n")
    cat("  sample_prefix: Prefix for output files\n")
    quit(status = 1)
  }
  
  paf_file <- args[1]
  output_dir <- args[2]
  sample_prefix <- args[3]
  
  # Check input file
  if (!file.exists(paf_file)) {
    cat("Error: PAF file not found:", paf_file, "\n")
    quit(status = 1)
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Parse PAF file
  alignments <- parse_paf(paf_file)
  
  if (nrow(alignments) == 0) {
    cat("No alignments found in PAF file\n")
    quit(status = 1)
  }
  
  # Get unique target contigs
  target_contigs <- unique(alignments$target_name)
  cat("Found", length(target_contigs), "target contigs\n")
  
  # Create single multi-page PDF with one page per target contig
  output_file <- file.path(output_dir, paste0(sample_prefix, "_alignment_plots.pdf"))
  cat("Creating multi-page PDF:", output_file, "\n")
  
  pdf(output_file, width = 12, height = 8)
  
  # Create plots for each target contig and add to PDF
  plots_created <- 0
  for (target_contig in target_contigs) {
    plot <- create_diagonal_plot(alignments, target_contig)
    if (!is.null(plot)) {
      print(plot)  # Print plot to PDF
      plots_created <- plots_created + 1
    }
  }
  
  dev.off()
  
  cat("Multi-page PDF created with", plots_created, "pages\n")
  cat("Saved:", output_file, "\n")
  cat("Alignment plotting completed successfully!\n")
}

# Run main function
if (!interactive()) {
  main()
} 
