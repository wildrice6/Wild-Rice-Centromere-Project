# ----------------------------------------------------------------------
# 1. Environment Setup and Library Loading
# ----------------------------------------------------------------------
# Ensure ggplot2 is installed; if not, please run: install.packages("ggplot2")
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# --- 2. Color Palette and Label Mapping Definitions ---

# Chromosome color mapping (Colorblind-friendly palette, adjusted from the provided list)
chromosome_colors <- c(
  "Chr01" = "#1F78B4", # Blue
  "Chr02" = "#A6CEE3", # Light Blue
  "Chr03" = "#FFBB78", # Orange
  "Chr04" = "#98DF8A", # Light Green
  "Chr05" = "#FB9A99", # Pink
  "Chr06" = "#CAB2D6", # Purple
  "Chr07" = "#A65628", # Brown
  "Chr08" = "#F781BF", # Deep Pink
  "Chr09" = "#8C8C8C", # Gray
  "Chr10" = "#BCBD22", # Yellow-Green
  "Chr11" = "#4DC4C4", # Teal
  "Chr12" = "#9EDAE5"  # Light Teal
)

# X-axis label mapping (retained unchanged)
label_mapping <- c(
  "AA_Ogla" = "O. glaberrima",
  "AA_Oruf" = "O. rufipogon",
  "AA_Oniv" = "O. nivara",
  "AA_Olon" = "O. longistaminata",
  "AA_Oglu" = "O. glumaepatula",
  "BB_Opun" = "O. punctata",
  "CC_Ooff" = "O. officinalis",
  "EE_Oaus" = "O. australiensis",
  "FF_Obra" = "O. brachyantha",
  "GG_Omey" = "O. meyeriana"
)


# ----------------------------------------------------------------------
# 3. Data Processing and Plotting Function Definitions
# ----------------------------------------------------------------------

process_and_plot <- function(input_path, output_path) {
  
  # Check if the input file exists
  if (!file.exists(input_path)) {
    stop(paste("Error: Input file does not exist:", input_path))
  }
  
  # 1. Load data
  cat(paste("Reading file:", input_path, "\n"))
  data <- read.table(
    input_path, 
    header = FALSE, 
    sep = "\t", 
    stringsAsFactors = FALSE,
    col.names = c("SeqA", "SeqB", "Similarity")
  )
  
  # 2. Data processing and feature extraction
  
  # Extract species names (first two fields of SeqA)
  data <- data %>%
    mutate(
      Species_Raw = str_extract(SeqA, "^[^_]+_[^_]+"),
      # Extract chromosome number (last field of SeqA)
      Chromosome = str_extract(SeqA, "Chr[0-9XYM]+$") 
    )
  
  # Map species names to readable labels
  data$Species_Label <- factor(data$Species_Raw, 
                               levels = names(label_mapping),
                               labels = label_mapping)
  
  # 3. Sort X-axis factors to ensure the plotting order matches label_mapping
  # (Species not present in label_mapping will be excluded from the factor levels)
  data$Species_Label <- factor(data$Species_Label, 
                               levels = unique(data$Species_Label))
  
  # 4. Generate scatter plot
  cat("Generating scatter plot...\n")
  
  p <- ggplot(data, aes(x = Species_Label, y = Similarity, color = Chromosome)) +
    
    # Scatter plot layer
    geom_point(size = 3, alpha = 0.7) + 
    
    # Apply custom colors
    scale_color_manual(
      values = chromosome_colors,
      name = "Chromosome"
    ) +
    
    # Set theme and axis labels
    labs(
      title = "Inter-Haplotype K-mer Jaccard Similarity by Chromosome",
      x = "Species/Genotype",
      y = "Jaccard Similarity"
    ) +
    
    # Adjust X-axis labels to vertical orientation for readability
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      panel.grid.major.x = element_blank(), # Remove vertical grid lines
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
  
  # 5. Save plot
  # High-quality PDF or PNG formats are recommended
  output_extension <- tools::file_ext(output_path)
  
  if (output_extension %in% c("png", "jpg", "jpeg", "tiff")) {
    ggsave(output_path, plot = p, width = 10, height = 6, dpi = 300)
  } else if (output_extension == "pdf") {
    ggsave(output_path, plot = p, width = 10, height = 6)
  } else {
    # Default to PDF
    output_path <- paste0(output_path, ".pdf")
    ggsave(output_path, plot = p, width = 10, height = 6)
    cat(paste("Warning: Unrecognized output format, defaulting to PDF:", output_path, "\n"))
  }
  
  cat(paste("Plot successfully saved to:", output_path, "\n"))
}

# ----------------------------------------------------------------------
# 4. Command-Line Argument Processing
# ----------------------------------------------------------------------
# Standard argument processing in R typically relies on the CLI environment or specific packages (e.g., optparse).
# A simplified method for retrieving command-line arguments is utilized here.

args <- commandArgs(trailingOnly = TRUE)

# Ensure the correct number of arguments (requires --input and --output)
if (length(args) == 0 || (length(args) != 4 && length(args) != 2)) {
  # Check for full --input INPUT --output OUTPUT structure
  stop("Usage: Rscript plot_similarity.R --input <input_file.tsv> --output <output_file.pdf>")
}

# Parse arguments
input_file <- NULL
output_file <- NULL

if (length(args) >= 2) {
  for (i in 1:length(args)) {
    if (args[i] == "--input" && i + 1 <= length(args)) {
      input_file <- args[i + 1]
    }
    if (args[i] == "--output" && i + 1 <= length(args)) {
      output_file <- args[i + 1]
    }
  }
}

if (is.null(input_file) || is.null(output_file)) {
  stop("Please provide paths for --input and --output files.")
}

# Execute the main function
process_and_plot(input_file, output_file)