#!/usr/bin/env Rscript

# Plot scatter plot comparing centromere sizes between haplotypes
#
# Usage:
#   Rscript plot_haplotype_centromere_scatter.R --input centromere_sizes.tsv --output scatter.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Plot haplotype centromere size comparison")
parser$add_argument("--input", required = TRUE, help = "Input TSV file")
parser$add_argument("--output", required = TRUE, help = "Output plot file (PDF/PNG)")
parser$add_argument("--width", type = "double", default = 8, 
                    help = "Plot width in inches (default: 8)")
parser$add_argument("--height", type = "double", default = 8, 
                    help = "Plot height in inches (default: 8)")
parser$add_argument("--point_size", type = "double", default = 3, 
                    help = "Point size (default: 3)")

args <- parser$parse_args()

# Read data
cat("[INFO] Reading input file:", args$input, "\n")
data <- read.table(args$input, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                   col.names = c("Chromosome", "CEN_Size", "Assembly"))

cat("[INFO] Found", nrow(data), "records\n")

# Parse assembly name to extract species and haplotype
data <- data %>%
  mutate(
    Assembly_parts = strsplit(Assembly, "_"),
    Species_code = sapply(Assembly_parts, function(x) {
      if (length(x) >= 2) paste(x[1], x[2], sep = "_") else x[1]
    }),
    Haplotype = sapply(Assembly_parts, function(x) {
      if (length(x) >= 3) x[length(x)] else "unknown"
    }),
    Assembly_parts = NULL  # Remove temporary column
  )

# Define species information with colors
species_info <- data.frame(
  Species_code = c("AA_Ogla", "AA_Oruf", "AA_Oniv", "AA_Olon", "AA_Oglu",
                   "BB_Opun", "CC_Ooff", "EE_Oaus", "FF_Obra", "GG_Omey"),
  color = c("#76D273", "#215A20", "#3BA738", "#51C54E", "#3D8347",
            "#F2AE2C", "#684E94", "#4E84C3", "#D55F6F", "#9D5427"),
  species_name = c("O. glaberrima", "O. rufipogon", "O. nivara",
                   "O. longistaminata", "O. glumaepatula", "O. punctata",
                   "O. officinalis", "O. australiensis", "O. brachyantha",
                   "O. meyeriana"),
  stringsAsFactors = FALSE
)

# Define species order for legend (top to bottom)
species_order <- c(
  "O. glaberrima",
  "O. rufipogon",
  "O. nivara",
  "O. longistaminata",
  "O. glumaepatula",
  "O. punctata",
  "O. officinalis",
  "O. australiensis",
  "O. brachyantha",
  "O. meyeriana"
)

# Define chromosome shapes mapping
chromosome_shapes_final <- c(
  "Chr01" = 1,
  "Chr02" = 2,
  "Chr03" = 3,
  "Chr04" = 4,
  "Chr05" = 5,
  "Chr06" = 6,
  "Chr07" = 7,
  "Chr08" = 8,
  "Chr09" = 9,
  "Chr10" = 10,
  "Chr11" = 11,
  "Chr12" = 12
)

# Add species information
data <- data %>%
  left_join(species_info, by = "Species_code") %>%
  mutate(species_name = ifelse(is.na(species_name), Species_code, species_name))

cat("[INFO] Found", length(unique(data$Species_code)), "species and", 
    length(unique(data$Haplotype)), "haplotypes\n")

# Check if we have exactly 2 haplotypes
haplotypes <- unique(data$Haplotype)
if (length(haplotypes) != 2) {
  cat("[WARN] Expected 2 haplotypes, found:", length(haplotypes), "\n")
  cat("      Haplotypes:", paste(haplotypes, collapse = ", "), "\n")
}

# Reshape data to wide format (one row per chromosome with hap1 and hap2 columns)
data_wide <- data %>%
  select(Chromosome, Species_code, species_name, color, Haplotype, CEN_Size) %>%
  pivot_wider(
    names_from = Haplotype,
    values_from = CEN_Size,
    names_prefix = "hap"
  )

# Check which haplotype columns exist
hap_cols <- grep("^hap", names(data_wide), value = TRUE)
cat("[INFO] Haplotype columns found:", paste(hap_cols, collapse = ", "), "\n")

# Use first two haplotype columns if they exist
if (length(hap_cols) < 2) {
  stop("[ERROR] Need at least 2 haplotypes for comparison")
}

hap1_col <- hap_cols[1]
hap2_col <- hap_cols[2]

# Rename for easier plotting
data_wide <- data_wide %>%
  rename(
    Hap1 = !!sym(hap1_col),
    Hap2 = !!sym(hap2_col)
  ) %>%
  filter(!is.na(Hap1) & !is.na(Hap2))

cat("[INFO] Comparing", hap1_col, "vs", hap2_col, "\n")
cat("[INFO]", nrow(data_wide), "chromosomes with both haplotypes\n")

# Convert species_name to factor with specified order
data_wide$species_name <- factor(data_wide$species_name, levels = species_order)

# Create named color vector (ordered by factor levels)
species_colors <- setNames(data_wide$color, data_wide$species_name)
species_colors <- species_colors[!duplicated(names(species_colors))]
# Reorder colors to match factor levels
species_colors <- species_colors[levels(data_wide$species_name)]
species_colors <- species_colors[!is.na(names(species_colors))]

# Convert to Mb for better readability
data_wide <- data_wide %>%
  mutate(
    Hap1_Mb = Hap1 / 1e6,
    Hap2_Mb = Hap2 / 1e6
  )

# Calculate axis limits (same for both axes)
max_val <- max(c(data_wide$Hap1_Mb, data_wide$Hap2_Mb), na.rm = TRUE)
min_val <- min(c(data_wide$Hap1_Mb, data_wide$Hap2_Mb), na.rm = TRUE)
axis_range <- max_val - min_val
axis_min <- max(0, min_val - 0.05 * axis_range)
axis_max <- max_val + 0.05 * axis_range

cat("[INFO] Creating scatter plot...\n")

# Create the plot with shape mapping
p <- ggplot(data_wide, aes(x = Hap1_Mb, y = Hap2_Mb, 
                           color = species_name, shape = Chromosome)) +
  geom_point(size = args$point_size, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_shape_manual(values = chromosome_shapes_final, name = "Chromosome") +
  scale_x_continuous(limits = c(axis_min, axis_max)) +
  scale_y_continuous(limits = c(axis_min, axis_max)) +
  labs(
    x = paste0("Centromere Size - ", hap1_col, " (Mb)"),
    y = paste0("Centromere Size - ", hap2_col, " (Mb)"),
    title = "Centromere Size Comparison Between Haplotypes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "italic"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    aspect.ratio = 1,
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4)),
    shape = guide_legend(order = 2, override.aes = list(size = 3))
  )

# Determine output format from file extension
output_ext <- tolower(tools::file_ext(args$output))

# Save plot
cat("[INFO] Saving plot to:", args$output, "\n")
if (output_ext == "pdf") {
  ggsave(args$output, plot = p, width = args$width, height = args$height, 
         device = "pdf")
} else if (output_ext %in% c("png", "jpg", "jpeg", "tiff")) {
  ggsave(args$output, plot = p, width = args$width, height = args$height, 
         dpi = 300)
} else {
  # Default to PDF
  cat("[WARN] Unknown file extension, saving as PDF\n")
  ggsave(args$output, plot = p, width = args$width, height = args$height, 
         device = "pdf")
}

# Print summary statistics
cat("\n[INFO] Summary statistics:\n")
summary_stats <- data_wide %>%
  group_by(species_name) %>%
  summarise(
    n_chromosomes = n(),
    mean_hap1 = mean(Hap1_Mb, na.rm = TRUE),
    mean_hap2 = mean(Hap2_Mb, na.rm = TRUE),
    cor_value = cor(Hap1_Mb, Hap2_Mb, use = "complete.obs"),
    .groups = "drop"
  )

for (i in 1:nrow(summary_stats)) {
  cat(sprintf("%s: n=%d, mean_%s=%.2f Mb, mean_%s=%.2f Mb, cor=%.3f\n",
              summary_stats$species_name[i],
              summary_stats$n_chromosomes[i],
              hap1_col, summary_stats$mean_hap1[i],
              hap2_col, summary_stats$mean_hap2[i],
              summary_stats$cor_value[i]))
}

# Overall correlation
overall_cor <- cor(data_wide$Hap1_Mb, data_wide$Hap2_Mb, use = "complete.obs")
cat(sprintf("\nOverall correlation: %.3f\n", overall_cor))

cat("\n[INFO] Done!\n")
