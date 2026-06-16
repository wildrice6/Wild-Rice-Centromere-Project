#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# 1. Load libraries
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
parser <- ArgumentParser(description = "Plot a grid of pie charts for chromosomal assembly components")
parser$add_argument("--input", required = TRUE, help = "Input TSV path")
parser$add_argument("--output", required = TRUE, help = "Output image path")
args <- parser$parse_args()

# -----------------------------------------------------------------------------
# 3. Define color palette
# -----------------------------------------------------------------------------
# Note: otherTE has been removed as it is now merged into 'others'
color_map <- c(
    'satellite'   = '#d62728',     
    'intactLTR'   = '#1f77b4',     
    'NUMT'        = '#ffbb78',     
    'rDNA'        = '#2ca02c',     
    'NUPT'        = '#8c564b',     
    'gene'        = '#9467bd',     
    'others'      = '#c7c7c7'      
)

# -----------------------------------------------------------------------------
# 4. Data reading and processing
# -----------------------------------------------------------------------------
cat("Processing data...\n")
raw_data <- read_tsv(args$input, col_types = cols())

# A. Data cleaning: processing percentages and restructuring categories
data_cleaned <- raw_data %>%
  mutate(
    pct_val = as.numeric(str_remove(Percentage, "%")),
    # [Modification]: Categorize all elements starting with "other_" as "others"
    ElementType = case_when(
      str_detect(ElementType, "^other_") ~ "others",
      TRUE ~ ElementType
    )
  )

# B. Obtain total value for 'Allelement' (representing total proportion of known functional elements)
all_element_totals <- data_cleaned %>%
  filter(ElementType == "Allelement") %>%
  select(assembly, Chromosome, total_pct = pct_val)

# C. Aggregate major components (excluding Allelement itself)
# Note: At this stage, "others" may already include sequences originally starting with "other_"
plot_data_main <- data_cleaned %>%
  filter(ElementType != "Allelement") %>%
  group_by(assembly, Chromosome, ElementType) %>%
  summarise(pct_val = sum(pct_val), .groups = "drop")

# D. Calculate "background" others (i.e., 100% minus Allelement)
# Then merge it with the existing "others" (those derived from other_) in plot_data_main
final_data <- plot_data_main %>%
  full_join(all_element_totals, by = c("assembly", "Chromosome")) %>%
  mutate(
    # If the type is others, its value is (accumulated other_ value) + (100 - total Allelement value)
    # If the type is not others, retain the original value
    pct_val = if_else(
      ElementType == "others",
      pct_val + pmax(0, 100 - total_pct),
      pct_val
    )
  ) %>%
  # Supplementary processing: if a row originally lacked other_ sequences, full_join may result in missing others rows
  # We must ensure background others are always present
  group_by(assembly, Chromosome) %>%
  group_modify(~ {
    if (!"others" %in% .x$ElementType) {
      bg_others <- pmax(0, 100 - .x$total_pct[1])
      return(bind_rows(.x, tibble(ElementType = "others", pct_val = bg_others)))
    }
    return(.x)
  }) %>%
  ungroup() %>%
  select(assembly, Chromosome, ElementType, pct_val)

# E. Sorting logic
# 1. Chromosome sorting
final_data <- final_data %>%
  mutate(chr_num = as.numeric(str_extract(Chromosome, "\\d+"))) %>%
  arrange(chr_num) %>%
  mutate(Chromosome = factor(Chromosome, levels = unique(Chromosome)))

# 2. Assembly sorting
final_data <- final_data %>%
  arrange(assembly) %>%
  mutate(assembly = factor(assembly, levels = unique(assembly)))

# 3. Element type sorting (following the order of color_map)
final_data <- final_data %>%
  mutate(ElementType = factor(ElementType, levels = names(color_map))) %>%
  filter(!is.na(ElementType)) # Filter out categories not defined in color_map

# -----------------------------------------------------------------------------
# 5. Plotting
# -----------------------------------------------------------------------------
cat("Generating pie chart grid...\n")

p <- ggplot(final_data, aes(x = "", y = pct_val, fill = ElementType)) +
  geom_bar(stat = "identity", position = "fill", width = 1, color = "white", linewidth = 0.05) +
  coord_polar(theta = "y", start = 0) +
  facet_grid(assembly ~ Chromosome, switch = "y") +
  scale_fill_manual(values = color_map, name = "Element Type") +
  labs(
    title = expression(paste("Centromere sequence composition in ", italic("Oryza"))),
    x = NULL, y = NULL
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b=15)),
    strip.text.x = element_text(size = 9, margin = margin(t=5, b=5)),
    strip.text.y.left = element_text(size = 9, angle = 0, hjust = 1, face = "italic"),
    legend.position = "right",
    panel.spacing = unit(0.1, "lines"),
    plot.margin = margin(10, 10, 10, 10)
  )

# -----------------------------------------------------------------------------
# 6. Saving
# -----------------------------------------------------------------------------
n_asm <- n_distinct(final_data$assembly)
n_chr <- n_distinct(final_data$Chromosome)

ggsave(
  args$output, 
  plot = p, 
  width = n_chr * 1.2 + 2, 
  height = n_asm * 0.8 + 1, 
  dpi = 300, 
  bg = "white"
)

cat(paste("Successfully saved to:", args$output, "\n"))