## Clear the environment
rm(list=ls())

## Set working directory
setwd("/red/mresende/viannam/PHG/output")

# Load required libraries
library(ggplot2)
library(dplyr)
library(data.table)

# Define directories
input_dir <- "/red/mresende/viannam/PHG/output"   # Directory containing assoc files
output_dir <- "/red/mresende/viannam/PHG/output" # Directory to save processed files

# Get the list of assoc files
file_list <- list.files(input_dir, pattern = "\\.pheno.gemma.assoc_updated.txt$", full.names = TRUE)

# Loop through each file
for (file_path in file_list) {
  # Extract base name without extension and "_updated"
  file_name <- basename(file_path)
  base_name <- sub("_GxE.pheno.gemma.assoc_updated\\.txt$", "", file_name)

  # Read data
  data <- fread(file_path)

  # Processar os dados diretamente com data.table
  gwas_data <- data[, .(
    marker_id = rs,
    chr = chr.x,
    #chr = as.integer(gsub("chr", "", chr.x, ignore.case = TRUE)), # Garantir que chr seja numérico
    ps,
    p_value = p_lrt,
    log_p = -log10(p_lrt)  # Calcular -log10(p-value)
  )][!is.na(chr)][order(chr, ps)]  # Filtrar NA e ordenar por chr e posição

  # Create Manhattan plot
  plot_title <- paste(base_name)
  threshold <- -log10(0.05/nrow(gwas_data))
  fdr_threshold <- -log10(0.05)
  output_file <- file.path(output_dir, paste0(base_name, ".jpg"))

  gwas_plot <- ggplot(gwas_data, aes(x = ps, y = log_p, color = as.factor(chr))) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = threshold, color = "black", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = fdr_threshold, color = "gray", linetype = "dotted", linewidth = 1) + # FDR threshold
    facet_wrap(~chr, scales = "free_x", nrow = 1, strip.position = "bottom") +  # Place chromosome labels on the x-axis
    theme_minimal() +
    labs(
      title = base_name,
      x = "Chromosome",
      y = "-log10(p-value)",
      color = "Chromosome"
    ) +
    scale_color_manual(values = rep(c("#1F77B4", "#4D4D4D"), length.out = length(unique(gwas_data$chr))), guide = "none") + # Custom alternating colors
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.background = element_blank(),  # Remove the background from facet labels
      strip.placement = "outside",         # Ensure facet labels are outside the plot area
      axis.text.x = element_blank(),       # Remove x-axis numbers
      axis.ticks.x = element_blank(),      # Remove x-axis ticks
      panel.spacing.x = unit(0.5, "lines") # Adjust spacing between panels
    )
  # Save the plot as a JPG file
  ggsave(output_file, plot = gwas_plot, width = 10, height = 6, dpi = 300)

  # Print message
  cat("Saved plot to:", output_file, "\n")
}

##-------------------retriving haplotypes------------------------------------##
# top_chr5 <- gwas_data %>%
#   filter(chr == 5) %>%       # Filter for chromosome 5
#   arrange(p_value) %>%       # Sort by p-value in ascending order
#   slice_head(n = 10)         # Select top 10 rows
#
# print(top_chr5)
