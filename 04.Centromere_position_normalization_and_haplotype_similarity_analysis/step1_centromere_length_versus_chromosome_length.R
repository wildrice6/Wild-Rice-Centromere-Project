#!/usr/bin/env Rscript

# 1. Load necessary libraries
if (!require("optparse", quietly = TRUE)) install.packages("optparse")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")

library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)

# 2. Argument parsing
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Path to input data file (TSV)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="facets_correlation_by_hap.pdf", help="Path to output PDF", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input file --input must be provided", call.=FALSE)
}

# 3. Read and process data
df <- read.table(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Species name mapping table
name_map <- c(
  "AA_Ogla" = "O. glaberrima", "AA_Oind" = "O. sativa ssp. indica", 
  "AA_Ojap" = "O. sativa ssp. japonica", "AA_Oruf" = "O. rufipogon", 
  "AA_Oniv" = "O. nivara", "AA_Olon" = "O. longistaminata", 
  "AA_Oglu" = "O. glumaepatula", "BB_Opun" = "O. punctata", 
  "CC_Ooff" = "O. officinalis", "EE_Oaus" = "O. australiensis", 
  "FF_Obra" = "O. brachyantha", "GG_Omey" = "O. meyeriana", 
  "XX_Lhex" = "Leersia hexandra"
)

# Specified species sorting order
ordered_species <- c(
  "O. sativa ssp. japonica", "O. sativa ssp. indica", "O. glaberrima", "O. rufipogon",
  "O. nivara", "O. longistaminata", "O. glumaepatula", "O. punctata",
  "O. officinalis", "O. australiensis", "O. brachyantha", "O. meyeriana"
)

# [Core Logic]: Construct facet labels and their ordering rules
# Desired order: Species 1(hap1), Species 1(hap2), Species 2(hap1), Species 2(hap2)...
facet_levels <- expand.grid(h = c("hap1", "hap2"), s = ordered_species) %>%
  mutate(full_label = paste0(s, "\n(", h, ")")) %>%
  pull(full_label)

# Data reshaping: Merge chr rows and cen rows
df_chr <- df %>% filter(type == "chr") %>% select(chrom, species, haplotype, length_chr = length)
df_cen <- df %>% filter(type == "cen") %>% select(chrom, species, haplotype, length_cen = length)

merged_df <- inner_join(df_chr, df_cen, by = c("chrom", "species", "haplotype")) %>%
  mutate(
    X_plot = length_chr / 1e6, # Mb
    Y_plot = length_cen / 1e3, # kb
    s_name = name_map[species],
    # Create facet labels: Species name + newline + haplotype
    facet_label = factor(paste0(s_name, "\n(", haplotype, ")"), levels = facet_levels)
  ) %>%
  # Filter out species not in the mapping list
  filter(!is.na(facet_label))

# 4. Calculate statistics for each facet (species + haplotype)
stats_df <- merged_df %>%
  group_by(facet_label) %>%
  summarize(
    r_val = cor(X_plot, Y_plot, use="complete.obs"),
    p_val = cor.test(X_plot, Y_plot)$p.value,
    x_pos = max(X_plot, na.rm=TRUE),
    y_pos = max(Y_plot, na.rm=TRUE)
  ) %>%
  mutate(
    label = paste0("r = ", sprintf("%.3f", r_val), "\n", "p = ", sprintf("%.3f", p_val))
  )

# 5. Plotting
p <- ggplot(merged_df, aes(x = X_plot, y = Y_plot)) +
  # Regression line
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.5) +
  # Scatter points
  geom_point(color = "grey40", size = 1.5, alpha = 0.7) +
  # Statistical annotation
  geom_text(data = stats_df, aes(x = x_pos, y = y_pos, label = label), 
            hjust = 1, vjust = 1.1, size = 2.5, fontface = "italic") +
  # [Faceting]: Facet by created facet_label, with 4 columns
  facet_wrap(~facet_label, scales = "free", ncol = 4) +
  labs(x = "Chromosome length (Mb)", y = "Centromere length (kb)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    # Set strip.text to italic; note that fine-tuning might be needed due to the addition of (hap)
    strip.text = element_text(face = "italic", size = 9),
    axis.text = element_text(color = "black", size = 7),
    axis.title = element_text(size = 10),
    legend.position = "none"
  )

# 6. Save output
# Increase height to accommodate multiple subplots (up to 12 species * 2 = 24 plots)
ggsave(opt$output, p, width = 10, height = 12, device = "pdf")

cat(paste0("Processing complete! Faceted by haplotype. Output file: ", opt$output, "\n"))