# -----------------------------------------------------------------------------
# 1. Load Required Libraries
# -----------------------------------------------------------------------------
# argparse: For parsing command-line arguments
# ggplot2: For data visualization
# readr: For efficient data reading
library(argparse)
library(ggplot2)
library(readr)

# -----------------------------------------------------------------------------
# 2. Configure Command-line Argument Parsing
# -----------------------------------------------------------------------------
# Create a parser instance
parser <- ArgumentParser(description = "Plot block size distribution from a TSV file as a line graph")

# Define the --input argument
parser$add_argument("--input", 
                    required = TRUE, 
                    help = "Path to the input TSV file")

# Define the --output argument
parser$add_argument("--output", 
                    required = TRUE, 
                    help = "Path to the output image file (e.g., my_plot.png)")

# Parse arguments provided via the command line
args <- parser$parse_args()

# -----------------------------------------------------------------------------
# 3. Define Mapping Between Material Names and Colors
# -----------------------------------------------------------------------------
# Create a named vector for consistent color palette configuration
color_palette <- c(
  "AA_Osat_jap" = "#59AC6E",
  "AA_Osat_ind" = "#CBE54E",
  "AA_Ogla" = "#76D273",
  "AA_Oruf" = "#215A20",
  "AA_Oniv" = "#3BA738",
  "AA_Olon" = "#51C54E",
  "AA_Oglu" = "#3D8347",
  "BB_Opun" = "#F2AE2C",
  "CC_Ooff" = "#684E94",
  "EE_Oaus" = "#4E84C3",
  "FF_Obra" = "#D55F6F",
  "GG_Omey" = "#9D5427"
)

# -----------------------------------------------------------------------------
# 4. Data Import and Preprocessing
# -----------------------------------------------------------------------------
# Retrieve input file path from arguments
input_file <- args$input

# Verify the existence of the input file
if (!file.exists(input_file)) {
  stop(paste("Error: Input file does not exist ->", input_file))
}

# Read data using read_table
col_names <- c("material_name", "block_size", "count", "percentage")
data <- read_table(input_file, col_names = col_names)

# Convert material_name to factor for proper categorical handling
data$material_name <- as.factor(data$material_name)

cat("Data successfully loaded from:", input_file, "\n")

# -----------------------------------------------------------------------------
# 5. Generate Line Plot Using ggplot2
# -----------------------------------------------------------------------------
cat("Generating visualization...\n")
p <- ggplot(data, aes(x = block_size, y = percentage, group = material_name, color = material_name)) +
  geom_line(size = 0.5) +
  # geom_point(size = 2.5) +
  scale_color_manual(values = color_palette) +
  labs(
    title = "Block Size Distribution by Material",
    x = "Block Size",
    y = "Percentage (%)",
    color = "Material Name"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

# -----------------------------------------------------------------------------
# 6. Export the Plot to a File
# -----------------------------------------------------------------------------
# Retrieve output file path from arguments
output_file <- args$output

# Save high-resolution image using ggsave
ggsave(output_file, plot = p, width = 10, height = 6, dpi = 300)

cat(paste("Plot successfully exported to:", output_file, "\n"))