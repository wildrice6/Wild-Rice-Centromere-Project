#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(argparse)
library(tidyr)
library(broom)

# Initialize argument parser
parser <- ArgumentParser(description='Plot chromosome similarity with permutation test significance (BH corrected)')
parser$add_argument('--input', required=TRUE, help='Input TSV file path')
parser$add_argument('--output', required=TRUE, help='Output plot file path')
parser$add_argument('--width', type='double', default=10, help='Plot width in inches (default 10)')
parser$add_argument('--height', type='double', default=8, help='Plot height in inches (default 8)')
parser$add_argument('--dpi', type='integer', default=300, help='Plot resolution (default 300)')
parser$add_argument('--perms', type='integer', default=10000, help='Number of permutations (default 10000)')
parser$add_argument('--sig_output_tsv', required=TRUE, help='Output TSV file path for significance values')

args <- parser$parse_args()

# --- Permutation Test Function ---
run_permutation_test_within_between <- function(df, n_perm = 10000) {
  within_vals <- df$avg_similarity[df$type == "within"]
  between_vals <- df$avg_similarity[df$type == "between"]

  if (length(within_vals) == 0 || length(between_vals) == 0) {
    return(NA)
  }

  obs_diff <- mean(within_vals) - mean(between_vals)
  all_vals <- c(within_vals, between_vals)
  n_within <- length(within_vals)
  
  # Perform permutations using replicate
  perm_diffs <- replicate(n_perm, {
    shuffled_vals <- sample(all_vals)
    perm_within <- shuffled_vals[1:n_within]
    perm_between <- shuffled_vals[(n_within + 1):length(all_vals)]
    mean(perm_within) - mean(perm_between)
  })

  # Calculate two-sided P-value
  # Note: High precision of raw P-values is essential for multiple testing correction
  p_val <- sum(abs(perm_diffs) >= abs(obs_diff)) / n_perm
  return(p_val)
}


# 1. Load data
cat("Reading file:", args$input, "\n")
raw_data <- read.table(args$input, sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names = NULL)

# 2. Extract species identifiers
raw_data <- raw_data %>%
  mutate(assembly_id = assembly) %>%
  mutate(species = sub("_hap[0-9]+$", "", assembly_id))

# 3. Compute significance and apply Benjamini-Hochberg (BH) correction
cat("Running permutation tests (n=", args$perms, ") per assembly...\n", sep="")

# Step A: Calculate raw P-values
significance_df <- raw_data %>%
  group_by(assembly_id) %>%
  do(data.frame(p_value = run_permutation_test_within_between(., n_perm = args$perms))) %>%
  ungroup() # Ungroup to apply BH correction across the entire column

# Step B: Perform BH correction
cat("Applying Benjamini-Hochberg correction...\n")
significance_df <- significance_df %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
  mutate(
    significance = case_when(
      adj_p_value < 0.001 ~ "FDR < 0.001", 
      adj_p_value < 0.01  ~ "FDR < 0.01",
      adj_p_value < 0.05  ~ "FDR < 0.05",
      is.na(adj_p_value)  ~ "Insufficient Data",
      TRUE                ~ "FDR >= 0.05"
    ),
    significance = factor(significance, levels = c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05", "FDR >= 0.05", "Insufficient Data"))
  )

# 4. Data processing: Aggregating means for visualization
cat("Processing data and calculating means...\n")
plot_data <- raw_data %>%
  group_by(assembly_id, species) %>%
  summarise(
    intra_mean = mean(avg_similarity[type == "within"], na.rm = TRUE),
    inter_mean = mean(avg_similarity[type == "between"], na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    intra_mean = intra_mean * 100,
    inter_mean = inter_mean * 100
  ) %>%
  filter(!is.na(intra_mean) & !is.na(inter_mean))

# 5. Merge significance metadata
plot_data <- plot_data %>%
  left_join(significance_df, by = "assembly_id")

# Save significance results (containing raw P-values and adjusted FDR)
cat("Saving significance data to TSV file:", args$sig_output_tsv, "\n")
sig_output_data <- plot_data %>%
  select(assembly_id, species, intra_mean, inter_mean, p_value, adj_p_value, significance)
write.table(sig_output_data, file = args$sig_output_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# 6. Define color scheme for species (consistent with established standards)
species_info <- data.frame(
  species = c("AA_Ojap", "AA_Oind", "AA_Ogla", "AA_Oruf", "AA_Oniv", "AA_Olon",
              "AA_Oglu", "BB_Opun", "CC_Ooff", "EE_Oaus", "FF_Obra", "GG_Omey", "XX_Lhex"),
  color = c("#59AC6E", "#CBE54E", "#76D273", "#215A20", "#3BA738", "#51C54E",
            "#3D8347", "#F2AE2C", "#684E94", "#4E84C3", "#D55F6F", "#9D5427", "#595959"),
  label = c("O. sativa ssp. japonica", "O. sativa ssp. indica", "O. glaberrima",
            "O. rufipogon", "O. nivara", "O. longistaminata", "O. glumaepatula",
            "O. punctata", "O. officinalis", "O. australiensis", "O. brachyantha",
            "O. meyeriana", "L. hexandra"),
  stringsAsFactors = FALSE
)

plot_data <- plot_data %>%
  left_join(species_info, by = "species") %>%
  mutate(species_label = factor(label, levels = species_info$label))

color_vector <- setNames(species_info$color, species_info$label)

# 7. Calculate coordinate axis ranges
all_vals <- c(plot_data$inter_mean, plot_data$intra_mean)
axis_min <- max(0, min(all_vals, na.rm=TRUE) - 5)
axis_max <- min(100, max(all_vals, na.rm=TRUE) + 5)

# 8. Generate plot
cat("\nCreating scatter plot...\n")
# Define shape mapping based on FDR significance categories
shape_values <- c(
  "FDR < 0.001" = 16,  
  "FDR < 0.01"  = 15,  
  "FDR < 0.05"  = 17,  
  "FDR >= 0.05" = 1, 
  "Insufficient Data" = 4 
)

p <- ggplot(plot_data, aes(x=inter_mean, y=intra_mean,
                           color=species_label,
                           shape=significance)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="gray50", alpha=0.5) +
  geom_point(alpha=0.8, size=6, stroke=1.2) +
  scale_color_manual(values=color_vector, name="Species") +
  scale_shape_manual(name = "Significance (BH Corrected)", values = shape_values) +
  labs(
    title = "Chromosome Similarity Across Species",
    subtitle = paste0("Significance: Permutation test with BH correction (n=", args$perms, ")"),
    x = "Inter-chromosome Similarity (%)",
    y = "Intra-chromosome Similarity (%)"
  ) +
  coord_fixed(ratio=1, xlim=c(axis_min, axis_max), ylim=c(axis_min, axis_max)) +
  theme_bw(base_size=14) +
  theme(
    plot.title = element_text(hjust=0.5, face="bold", size=16),
    plot.subtitle = element_text(hjust=0.5, size=10),
    legend.position = "right",
    legend.title = element_text(face="bold")
  ) +
  guides(
    color = guide_legend(order=1, override.aes=list(size=4)),
    shape = guide_legend(order=2, override.aes = list(size = 4, stroke = 1.2, fill = NA))
  )

# Save the final plot
cat("Saving plot to:", args$output, "\n")
ggsave(args$output, plot=p, width=args$width, height=args$height, dpi=args$dpi)

cat("Done!\n")