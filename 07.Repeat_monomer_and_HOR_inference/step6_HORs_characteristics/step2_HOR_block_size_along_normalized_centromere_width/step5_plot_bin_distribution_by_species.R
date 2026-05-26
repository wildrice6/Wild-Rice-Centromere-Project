#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

# --- 1. Define Command Line Arguments ---
parser <- ArgumentParser(description = "Visualize block size distribution across genomic bins for various assemblies.")
parser$add_argument("-i", "--input", required = TRUE, help = "Path to the input TSV file.")
parser$add_argument("-o", "--output", required = TRUE, help = "Path for the output plot (e.g., plot.pdf or plot.png).")

args <- parser$parse_args()

# --- 2. Define Custom Color Palette ---
# A named vector mapping species names to specific hex color codes for consistency.
color_palette <- c(
  "O. sativa ssp. japonica" = "#59AC6E",
  "O. sativa ssp. indica" = "#CBE54E",
  "O. glaberrima" = "#76D273",
  "O. rufipogon" = "#215A20",
  "O. nivara" = "#3BA738",
  "O. longistaminata" = "#51C54E",
  "O. glumaepatula" = "#3D8347",
  "O. punctata" = "#F2AE2C",
  "O. officinalis" = "#684E94",
  "O. australiensis" = "#4E84C3",
  "O. brachyantha" = "#D55F6F",
  "O. meyeriana" = "#9D5427"
)

cat("--- Reading input file:", args$input, "---\n")

# --- 3. Data Loading and Preprocessing ---
plot_data <- read_tsv(args$input, col_types = cols())

# Generate a 'species' column by matching the 'assembly' column 
# against the species identifiers defined in the color palette.
plot_data <- plot_data %>%
  mutate(species = case_when(
    # str_detect identifies if the assembly name string contains the species identifier
    str_detect(assembly, "O. sativa ssp. japonica") ~ "O. sativa ssp. japonica",
    str_detect(assembly, "O. sativa ssp. indica") ~ "O. sativa ssp. indica",
    str_detect(assembly, "O. glaberrima") ~ "O. glaberrima",
    str_detect(assembly, "O. rufipogon") ~ "O. rufipogon",
    str_detect(assembly, "O. nivara") ~ "O. nivara",
    str_detect(assembly, "O. longistaminata") ~ "O. longistaminata",
    str_detect(assembly, "O. glumaepatula") ~ "O. glumaepatula",
    str_detect(assembly, "O. punctata") ~ "O. punctata",
    str_detect(assembly, "O. officinalis") ~ "O. officinalis",
    str_detect(assembly, "O. australiensis") ~ "O. australiensis",
    str_detect(assembly, "O. brachyantha") ~ "O. brachyantha",
    str_detect(assembly, "O. meyeriana") ~ "O. meyeriana",
    TRUE ~ as.character(assembly) # Fallback to original name if no match is found
  ))

# Identify and report any assemblies that failed to match the predefined palette
unmatched_species <- setdiff(plot_data$species, names(color_palette))
if (length(unmatched_species) > 0) {
  cat("Warning: The following assemblies did not match the species palette and will use default colors:\n")
  cat(paste(unmatched_species, collapse = "\n"), "\n")
}

cat("--- Data processing complete. Generating plot... ---\n")

# --- 4. Data Visualization using ggplot2 ---
p <- ggplot(plot_data, aes(x = bin_number, y = mean_average_max_block_size, group = assembly, color = species)) +
  # Lines are grouped by assembly to ensure distinct lines for individual haplotypes/assemblies
  # while sharing the same color based on the parent species.
  geom_line(linewidth = 0.8, alpha = 0.9) +
  
  # Apply the predefined color mapping
  scale_color_manual(values = color_palette) +
  
  # Title and axis labeling
  labs(
    title = "Mean Average Max Block Size Across Bins",
    subtitle = "Lines categorized by species; individual lines represent distinct assemblies",
    x = "Bin Number",
    y = "Mean Average Max Block Size",
    color = "Species"
  ) +
  
  # Apply a minimal clean theme
  theme_classic() +
  
  # Refine theme parameters for publication-quality output
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

# --- 5. Export Output ---
cat("--- Saving plot to:", args$output, "---\n")

# ggsave handles format selection based on file extension (.pdf, .png, .svg)
ggsave(
  filename = args$output, 
  plot = p, 
  width = 12,   # inches
  height = 7,   # inches
  dpi = 300     # resolution
)

cat("--- Execution completed successfully. ---\n")