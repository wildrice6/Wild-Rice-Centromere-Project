library(ggplot2)
library(dplyr)
library(viridis)
library(cowplot)

# ==========================================
# 1. Load LTR Homology Data
# ==========================================
ltr_file <- "D:/LTR_Species_Identity_Stats.csv"
ltr_data <- read.csv(ltr_file, stringsAsFactors = FALSE)

# Retain only required columns
ltr_data <- ltr_data %>%
  select(Species, Percent_gt_097) %>%
  rename(species_code = Species, identity_percent = Percent_gt_097)

cat("LTR homology data loaded successfully:\n")
print(ltr_data)

# ==========================================
# 2. Configure HOR Analysis Paths and Species Mapping
# ==========================================
base_path <- "D:/HOR_score_results/"
crm_path <- paste0(base_path, "CRM_intervals_HOR_scores.tsv")
rand_path <- paste0(base_path, "Random_intervals_HOR_scores.tsv")

output_pdf <- paste0(base_path, "Species_Level_HOR_with_LTR_Bubbles.pdf")

# Species name mapping table
name_map <- c(
  "AA_Osat_jap" = "O. sativa ssp. japonica",
  "AA_Osat_ind" = "O. sativa ssp. indica",
  "AA_Ogla"     = "O. glaberrima",
  "AA_Oruf"     = "O. rufipogon",
  "AA_Oniv"     = "O. nivara",
  "AA_Olon"     = "O. longistaminata",
  "AA_Oglu"     = "O. glumaepatula",
  "BB_Opun"     = "O. punctata",
  "CC_Ooff"     = "O. officinalis",
  "EE_Oaus"     = "O. australiensis"
)

# Defined display order
species_order <- unname(name_map)

# ==========================================
# 3. Process HOR Data
# ==========================================
prepare_data_species <- function(file_path, group_label) {
  df <- read.table(file_path, header=TRUE, sep="\t", check.names=FALSE, fill=TRUE, quote="")
  
  # Extract the first column (ID) and the last column (Score)
  df_clean <- df[, c(1, ncol(df))]
  colnames(df_clean) <- c("id", "raw_score")
  
  df_clean %>%
    mutate(sp_id = case_when(
      grepl("AA_Osat_jap", id) ~ "AA_Osat_jap",
      grepl("AA_Osat_ind", id) ~ "AA_Osat_ind",
      grepl("AA_Ogla", id) ~ "AA_Ogla",
      grepl("AA_Oruf", id) ~ "AA_Oruf",
      grepl("AA_Oniv", id) ~ "AA_Oniv",
      grepl("AA_Olon", id) ~ "AA_Olon",
      grepl("AA_Oglu", id) ~ "AA_Oglu",
      grepl("BB_Opun", id) ~ "BB_Opun",
      grepl("CC_Ooff", id) ~ "CC_Ooff",
      grepl("EE_Oaus", id) ~ "EE_Oaus",
      TRUE ~ "Other"
    )) %>%
    filter(sp_id != "Other") %>%
    mutate(species = name_map[sp_id]) %>%
    mutate(score = as.numeric(as.character(raw_score))) %>%
    mutate(score = ifelse(is.na(score), 0, score)) %>%
    select(species, score) %>%
    mutate(Group = group_label)
}

# Merge CRM and Random datasets
combined_df <- rbind(
  prepare_data_species(crm_path, "CRM"), 
  prepare_data_species(rand_path, "Random")
)

# Set factor levels for categorical variables
combined_df$species <- factor(combined_df$species, levels = species_order)

# ==========================================
# 4. Statistical Analysis and Summary Table
# ==========================================
stats_summary <- data.frame()
set.seed(123)

for (sp in species_order) {
  sub_data <- combined_df %>% filter(species == sp)
  crm_s <- sub_data$score[sub_data$Group == "CRM"]
  rnd_s <- sub_data$score[sub_data$Group == "Random"]
  
  if(length(crm_s) > 0 && length(rnd_s) > 0) {
    m_crm <- mean(crm_s); m_rnd <- mean(rnd_s)
    diff <- m_rnd - m_crm
    
    # Permutation test
    all_s <- sub_data$score
    n_crm <- length(crm_s)
    perm_diffs <- replicate(1000, {
      shuffled <- sample(all_s)
      mean(shuffled[(n_crm+1):length(shuffled)]) - mean(shuffled[1:n_crm])
    })
    p_val <- (sum(perm_diffs >= diff) + 1) / 1001
    sig <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else if(p_val < 0.05) "*" else "ns"
    
    stats_summary <- rbind(stats_summary, data.frame(
      Full_Name = sp,
      CRM_Mean = round(m_crm, 4),
      Random_Mean = round(m_rnd, 4),
      Difference = round(diff, 4),
      P_Value = round(p_val, 4),
      Significance = sig
    ))
  }
}

# ==========================================
# 5. Prepare LTR Bubble Data
# ==========================================
# Convert species codes in LTR data to full scientific names
ltr_data$species <- name_map[ltr_data$species_code]
ltr_data <- ltr_data[!is.na(ltr_data$species), ]

# Ensure species order consistency with the HOR plot
ltr_data$species <- factor(ltr_data$species, levels = species_order)

# Map identity percentages to an 80-100% range (cool to warm tones)
ltr_data <- ltr_data %>%
  mutate(
    # Normalize percentage to 0-1 range (based on 80-100% scale)
    color_value = (identity_percent - 80) / 20,
    # Constrain values within the 0-1 interval
    color_value = pmax(0, pmin(1, color_value))
  )

# Calculate bubble coordinates
max_hor_score <- max(combined_df$score)
star_y <- max_hor_score * 1.05
bubble_y <- max_hor_score * 1.15
y_upper_limit <- max_hor_score * 1.25

# Generate bubble plot dataset
bubble_plot_data <- ltr_data %>%
  mutate(
    x_pos = as.numeric(species),
    y_pos = bubble_y,
    # Bubble size: apply a narrow range (8-10) for subtle differentiation
    size = 8 + (identity_percent - min(identity_percent)) / 
      (max(identity_percent) - min(identity_percent)) * 2
  )

# ==========================================
# 6. Create Main Plot (Manual Annotation to Avoid Aesthetic Mapping Conflicts)
# ==========================================
# Initialize base plot (boxplots and jittered points)
base_plot <- ggplot() +
  # 1. HOR boxplots
  geom_boxplot(
    data = combined_df,
    aes(x = species, y = score, fill = Group),
    outlier.shape = NA,
    alpha = 0.6,
    width = 0.6,
    position = position_dodge(0.8),
    color = "black",
    size = 0.3
  ) +
  # 2. HOR jittered points
  geom_point(
    data = combined_df,
    aes(x = species, y = score, color = Group),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1,
    alpha = 0.5
  ) +
  # 3. Significance asterisks
  geom_text(
    data = stats_summary,
    aes(x = Full_Name, y = star_y, label = Significance),
    size = 5,
    fontface = "bold",
    vjust = 0.5
  )

# Manually annotate bubbles to prevent 'fill' aesthetic conflicts
bubble_colors <- viridis(100, option = "plasma")

for(i in 1:nrow(bubble_plot_data)) {
  row <- bubble_plot_data[i, ]
  
  # Retrieve color index
  color_idx <- round(row$color_value * 99) + 1
  color_idx <- min(max(color_idx, 1), 100)
  
  base_plot <- base_plot +
    annotate(
      "point",
      x = row$x_pos,
      y = row$y_pos,
      size = row$size,
      fill = bubble_colors[color_idx],
      color = "black",
      shape = 21,
      stroke = 0.8,
      alpha = 0.8
    )
}

# ==========================================
# 7. Add Scale, Labels, and Theme Styles
# ==========================================
main_plot <- base_plot +
  # Fill colors for HOR boxplots
  scale_fill_manual(
    name = "Region Type",
    values = c("CRM" = "#1F77B4", "Random" = "#FF7F0E"),
    guide = guide_legend(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 0.7)
    )
  ) +
  # Colors for HOR jittered points
  scale_color_manual(
    name = "Region Type",
    values = c("CRM" = "#1F77B4", "Random" = "#FF7F0E"),
    guide = "none"
  ) +
  # Y-axis range configuration
  scale_y_continuous(
    limits = c(0, y_upper_limit),
    expand = expansion(mult = c(0, 0.15))
  ) +
  # Labels and titles
  labs(
    title = "HOR Score Comparison with LTR Retrotransposon Identity",
    subtitle = "Box plots: HOR scores | Top bubbles: LTR identity ≥97% (cool:80% to warm:100%)",
    x = "Oryza Species",
    y = "HOR Score",
    caption = "P-values from 1000-fold permutation test: *** p < 0.001, ** p < 0.01, * p < 0.05"
  ) +
  # Theme configuration
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "bold.italic",
      size = 11
    ),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(
      size = 13,
      face = "bold",
      margin = margin(t = 15)
    ),
    axis.title.y = element_text(
      size = 13,
      face = "bold",
      margin = margin(r = 15)
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 10)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 12,
      margin = margin(b = 20)
    ),
    plot.caption = element_text(
      hjust = 0.5,
      size = 10,
      color = "gray40",
      margin = margin(t = 15)
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 11, face = "bold", margin = margin(b = 5)),
    legend.text = element_text(size = 10),
    legend.spacing.y = unit(0.4, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# ==========================================
# 8. Create Color Gradient Legend
# ==========================================
# Generate data for bubble color legend
color_legend_data <- data.frame(
  value = seq(0, 1, length.out = 100),
  x = 1:100,
  y = 1
)

# Create color gradient legend plot
color_legend_plot <- ggplot(color_legend_data, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = viridis(100, option = "plasma"),
    name = "LTR Identity ≥ 97% (%)",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("80%", "90%", "100%"),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3, "cm"),
      barheight = unit(0.4, "cm"),
      direction = "horizontal",
      title.theme = element_text(size = 11, face = "bold"),
      label.theme = element_text(size = 9)
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(5, 5, 5, 5)
  )

# ==========================================
# 9. Extract and Combine Legends
# ==========================================
# Extract the legend for Region Type
main_legend <- get_legend(main_plot)

# Extract the color gradient legend
color_legend <- get_legend(color_legend_plot)

# Create a blank plot for spacing
blank_plot <- ggplot() + theme_void()

# Construct the right-side legend panel (vertical arrangement)
right_legend_panel <- plot_grid(
  main_legend,
  blank_plot,
  color_legend,
  ncol = 1,
  rel_heights = c(0.4, 0.1, 0.4),
  align = "v"
)

# Remove original legend from the main plot
main_plot_no_legend <- main_plot + theme(legend.position = "none")

# ==========================================
# 10. Assemble Main Plot and Legend Panels
# ==========================================
final_plot <- plot_grid(
  main_plot_no_legend,
  right_legend_panel,
  ncol = 2,
  rel_widths = c(0.85, 0.15),
  align = "h"
)

# ==========================================
# 11. Export to PDF
# ==========================================
ggsave(
  output_pdf,
  plot = final_plot,
  width = 18,  # Increased width to accommodate legends
  height = 10,
  dpi = 300,
  device = "pdf"
)

cat("\n=== Analysis Complete ===\n")
cat("PDF file saved to:", output_pdf, "\n")
cat("\nVisualization Features:\n")
cat("1. HOR boxplots (Blue: CRM, Orange: Random)\n")
cat("2. Significance asterisks positioned above boxplots\n")
cat("3. Color-coded bubbles above asterisks (Cool: 80% to Warm: 100%)\n")
cat("4. Right panel includes legends for:\n")
cat("   • Region Type (boxplot colors)\n")
cat("   • LTR Identity Gradient (80%-100%)\n")
cat("\nBubble Configuration:\n")
cat("• Palette: viridis plasma\n")
cat("• Range: Cool tones (80%) → Warm tones (100%)\n")
cat("• Size: 8-10 (narrow variation)\n")
cat("• Numeric percentage labels are omitted\n")
cat("\nLTR Homology Statistics:\n")
print(ltr_data)

cat("\nHOR Statistical Summary:\n")
print(stats_summary)

# Render plot in the graphics device
print(final_plot)

--- END OF FILE ---