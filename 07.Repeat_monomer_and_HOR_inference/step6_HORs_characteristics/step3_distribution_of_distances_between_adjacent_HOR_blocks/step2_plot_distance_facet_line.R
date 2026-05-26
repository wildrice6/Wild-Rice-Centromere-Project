#!/usr/bin/env Rscript

# Load required packages
# If these packages are not yet installed, please run the following in the R console:
# install.packages(c("ggplot2", "readr", "argparse", "scales"))
library(ggplot2)
library(readr)
library(argparse)
library(scales)

# 1. Define command-line argument parser
parser <- ArgumentParser(description = "Generate faceted line plots from TSV files")
parser$add_argument("--input", type = "character", required = TRUE,
                    help = "Path to the input TSV file (file must contain a header)")
parser$add_argument("--output", type = "character", required = TRUE,
                    help = "Path to the output image file (e.g., my_plot.png)")

# Parse arguments
args <- parser$parse_args()

# 2. Define color palette
color_palette <- c(
    "AA_Osat_jap" = "#59AC6E",
    "AA_Osat_ind" = "#CBE54E",
    "AA_Ogla"     = "#76D273",
    "AA_Oruf"     = "#215A20",
    "AA_Oniv"     = "#3BA738",
    "AA_Olon"     = "#51C54E",
    "AA_Oglu"     = "#3D8347",
    "BB_Opun"     = "#F2AE2C",
    "CC_Ooff"     = "#684E94",
    "EE_Oaus"     = "#4E84C3",
    "FF_Obra"     = "#D55F6F",
    "GG_Omey"     = "#9D5427"
)

# 3. Data reading and processing
# --- Primary modifications ---
message(paste("Reading input file:", args$input))
# Modification: Column names are no longer manually specified; read_tsv automatically retrieves the header from the first row.
# read_tsv automatically infers data types for each column, correctly identifying the second column as numeric.
data <- read_tsv(args$input, show_col_types = FALSE)

# [Optional but recommended] For code robustness, columns can be checked or renamed
# to ensure consistency with the names used in subsequent ggplot code.
# Assuming the file headers are "material_name", "number_div_1M", "count", "percentage"
# If not, the following line can be used to enforce renaming:
# names(data) <- c("material_name", "number_div_1M", "count", "percentage")
# --------------------

# 4. Visualization using ggplot2
message("Generating plot...")
# Verify existence of required columns
if (!all(c("material_name", "number_div_1M", "percentage") %in% names(data))) {
  stop("Input file is missing required columns: 'material_name', 'number_div_1M', 'percentage'")
}

plot <- ggplot(data, aes(x = number_div_1M, y = percentage, color = material_name)) +
  geom_line(linewidth = 0.8) +
  
  facet_wrap(~ material_name, ncol = 4, scales = "free_y") +
  
  scale_color_manual(values = color_palette) +
  
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  scale_y_continuous(breaks = pretty_breaks(n = 4)) +
  
  labs(
    title = "Trend of Percentages Across Materials Relative to number_div_1M",
    x = "Number_div_1M",
    y = "Percentage (%)",
    color = "Material Name"
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "none", 
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

# 5. Save the plot
message(paste("Saving plot to:", args$output))
ggsave(
  filename = args$output,
  plot = plot,
  width = 12,
  height = 8,
  dpi = 300
)

message("Script execution complete!")