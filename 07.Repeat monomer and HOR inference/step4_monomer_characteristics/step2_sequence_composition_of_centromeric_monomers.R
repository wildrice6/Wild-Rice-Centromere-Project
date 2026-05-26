#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# 1. Load required libraries
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))     
suppressPackageStartupMessages(library(stringr))    
suppressPackageStartupMessages(library(ggplot2))

# -----------------------------------------------------------------------------
# 2. Argument parsing
# -----------------------------------------------------------------------------
parser <- ArgumentParser(description = "Plot a chromosome-assembly component pie chart grid from a TSV file")

parser$add_argument("--input", 
                    required = TRUE, 
                    help = "Path to the input TSV file")

parser$add_argument("--output", 
                    required = TRUE, 
                    help = "Path for the output image file (e.g., pie_grid.png)")

args <- parser$parse_args()

# -----------------------------------------------------------------------------
# 3. Define colors and read data
# -----------------------------------------------------------------------------
input_file <- args$input
output_file <- args$output

if (!file.exists(input_file)) {
  stop(paste("Error: Input file not found ->", input_file))
}

# [Modification: Removed color mapping for 'fragment']
color_map <- c(
    'satellite'   = '#d62728',     # Strong Red
    'intactLTR'   = '#1f77b4',     # Dark Blue
    'NUMT'        = '#ffbb78',     # Light Orange/Yellow
    'rDNA'        = '#2ca02c',     # Strong Green
    'NUPT'        = '#8c564b',     # Brown/Coffee
    'gene'        = '#9467bd',     # Purple
    'others'      = '#c7c7c7'      # Grey
)

cat("Reading data from", input_file, "...\n")
data <- read_tsv(input_file, col_types = cols())

# -----------------------------------------------------------------------------
# 4. Data preprocessing
# -----------------------------------------------------------------------------
cat("Processing data...\n")

# Define the expected order for species and haplotypes
species_base_order <- c(
  "O. sativa ssp. japonica", "O. sativa ssp. indica", "O. glaberrima",
  "O. rufipogon", "O. nivara", "O. longistaminata", "O. glumaepatula",
  "O. punctata", "O. officinalis", "O. australiensis", "O. brachyantha",
  "O. meyeriana", "L. hexandra"
)
hap_order <- c("hap1", "hap2", "hap3", "hap4")

# [Core Modification: Exclude 'chromosome', 'assembly', and 'fragment' when extracting element columns]
element_cols <- setdiff(names(data), c("chromosome", "assembly", "fragment"))

data_long <- data %>%
  # Ensure all element columns are numeric
  mutate(across(all_of(intersect(names(.), element_cols)), as.numeric)) %>%
  mutate(
    # Calculate the sum of known elements excluding 'fragment'
    known_sum = rowSums(select(., all_of(element_cols)), na.rm = TRUE),
    # Assign the remaining proportion (including the original 'fragment') to 'others'
    others = pmax(0, 1 - known_sum)
  ) %>%
  select(-known_sum) %>%
  # Exclude 'fragment' from the selection if it exists to avoid entering the long format table
  select(-any_of("fragment")) %>%
  pivot_longer(
    cols = c(all_of(element_cols), "others"), 
    names_to = "element_type", 
    values_to = "proportion"
  ) %>%
  filter(proportion > 0) %>%
  
  # 1. Use separate to split the 'assembly' column into species and haplotype parts
  separate(assembly, into = c("species_part", "hap_part"), sep = "_(?=[^_]+$)", remove = FALSE) %>%

  # 2. Convert the split parts into ordered factors
  mutate(
    species_part = factor(species_part, levels = species_base_order),
    hap_part = factor(hap_part, levels = hap_order)
  ) %>%
  
  # 3. Sort the entire data frame based on these two new ordered factors
  arrange(species_part, hap_part) %>%
  
  # 4. Apply natural sorting to the X-axis (chromosomes)
  mutate(chr_num = as.numeric(str_extract(chromosome, "\\d+"))) %>%
  arrange(chr_num) %>%
  
  # 5. Convert to factors to maintain display order after sorting
  mutate(
    assembly = factor(assembly, levels = unique(assembly)),
    chromosome = factor(chromosome, levels = unique(chromosome)),
    element_type = factor(element_type, levels = names(color_map))
  ) %>%
  
  # 6. Remove temporary columns
  select(-species_part, -hap_part, -chr_num)


cat("Data processing complete.\n")

# -----------------------------------------------------------------------------
# 5. Plot the pie chart grid using ggplot2
# -----------------------------------------------------------------------------
cat("Generating plot...\n")

p <- ggplot(data_long, 
            aes(x = "", y = proportion, fill = element_type)) +
  geom_bar(stat = "identity", position = "fill", width = 1, color = "white", linewidth = 0.2) +
  coord_polar(theta = "y", start = 0) +
  facet_grid(assembly ~ chromosome, drop = FALSE) +
  scale_fill_manual(values = color_map, name = "Element Type", drop = FALSE) +
  labs(title = "Genomic Component Proportions (Excl. Fragment)", x = NULL, y = NULL) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b=10)),
    strip.text.x = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# -----------------------------------------------------------------------------
# 6. Save plot to file
# -----------------------------------------------------------------------------
# Calculate dynamic height
n_assemblies <- n_distinct(data$assembly)
dynamic_height <- max(4, 1.0 * n_assemblies)

ggsave(output_file, plot = p, width = 16, height = dynamic_height, dpi = 300, bg = "white")

cat(paste("Plot successfully saved to:", output_file, "\n"))