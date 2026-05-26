#!/usr/bin/env Rscript

# Load required R packages
# If not installed, please run: install.packages(c("argparse", "readr", "dplyr", "ggplot2", "gtools"))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gtools))

# --- 1. Command-line Argument Parsing ---
parser <- ArgumentParser(description = "Generate violin and jitter plots of insertion times grouped by chromosome.")

parser$add_argument("--input",
                    required = TRUE,
                    help = "Path to the input tab- or space-separated data file.")

parser$add_argument("--output",
                    required = TRUE,
                    help = "Path to the output image file (e.g., plot.png, plot.pdf).")

# Parse arguments
args <- parser$parse_args()

# --- 2. Data Acquisition and Processing ---

# Define column names
column_names <- c(
  'chromosome', 'start', 'end', 'id', 
  'score', 'strand', 'insertion_time'
)

# Use readr::read_table to read data, which handles one or more spaces as delimiters
# Use tryCatch to handle potential errors such as missing or empty files
tryCatch({
  df <- readr::read_table(
    args$input, 
    col_names = column_names,
    col_types = cols(.default = "c", start="i", end="i", score="d", insertion_time="d") # Optimize reading efficiency
  )
  
  if (nrow(df) == 0) {
    stop("Input file is empty.")
  }

  # Data transformation: Convert insertion time to millions of years (mya)
  df <- df %>%
    mutate(insertion_time_mya = insertion_time / 1000000)
  
  # Perform natural sorting on chromosome names and set them as factor levels
  # This ensures ggplot2 renders the axes in the correct alphanumeric order
  sorted_chromosomes <- gtools::mixedsort(unique(df$chromosome))
  df$chromosome <- factor(df$chromosome, levels = sorted_chromosomes)

}, error = function(e) {
  cat(sprintf("Error: %s\n", e$message))
  quit(status = 1)
})


# --- 3. Color Definition and Visualization ---

# Define color palette for each chromosome (consistent with the Python implementation)
chromosome_colors <- c(
  "Chr01" = "#1F78B4",
  "Chr02" = "#A6CEE3",
  "Chr03" = "#FFBB78",
  "Chr04" = "#98DF8A",
  "Chr05" = "#FB9A99",
  "Chr06" = "#CAB2D6",
  "Chr07" = "#A65628",
  "Chr08" = "#F781BF",
  "Chr09" = "#8C8C8C",
  "Chr10" = "#BCBD22", 
  "Chr11" = "#4DC4C4", 
  "Chr12" = "#9EDAE5"  
)

# Construct visualization using ggplot2
p <- ggplot(df, aes(x = chromosome, y = insertion_time_mya)) +
  
  # 6. Plot violin diagrams with custom color mapping
  # aes(fill = chromosome) maps the fill aesthetic to the chromosome variable
  # trim=FALSE ensures the violin density is not truncated at the data extremes
  geom_violin(aes(fill = chromosome), trim = FALSE, show.legend = FALSE) +
  
  # 7. Overlay semi-transparent jittered points (geom_jitter)
  # shape = 21 specifies a fillable circle
  # fill = "grey" sets the point interior color
  # alpha = 0.1 controls transparency to manage overplotting
  geom_jitter(shape = 21, fill = "grey", size = 0.5, stroke = 0.7, alpha = 0.1, width=0.2) +
  
  # ==========================================================================
  # <<< Revision Point >>>
  # 8. Render three red dashed horizontal reference lines
  geom_hline(yintercept = 0.72, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 5.02, color = "red", linetype = "dashed", size = 1) + # <-- Newly added reference line
  geom_hline(yintercept = 15, color = "red", linetype = "dashed", size = 1) +
  # ==========================================================================
  
  # Apply the custom color palette
  scale_fill_manual(values = chromosome_colors) +
  
  # 9. Define plot title and axis labels
  labs(
    title = "Distribution of Insertion Time by Chromosome",
    x = "Chromosome",
    y = "Insertion Time (Million Years)"
  ) +
  
  # Apply a clean theme style similar to Seaborn's "whitegrid"
  theme_bw() +
  
  # Refine theme parameters
  theme(
    plot.title = element_text(size = 16, hjust = 0.5), # Center-align title
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # Rotate X-axis labels
    panel.grid.major.x = element_blank(), # Remove vertical major grid lines
    panel.grid.minor = element_blank()
  )

# --- 4. Export Figure ---
ggsave(
  args$output, 
  plot = p, 
  width = 14, 
  height = 8, 
  units = "in", 
  dpi = 300
)

cat(sprintf("Visualization complete! Image saved to '%s'\n", args$output))

--- END OF FILE ---