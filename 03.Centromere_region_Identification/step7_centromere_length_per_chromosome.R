# 1. Load Libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# 2. Set File Paths
file_path <- "C:/Users/L1.xlsx"

# 3. Data Cleaning and Preprocessing
df_clean <- read_excel(file_path) %>%
  mutate(Length = as.numeric(`Length (bp)`)) %>%
  # Remove all NA values
  filter(!is.na(Length)) %>%
  # Handle duplicate records to ensure a unique entry per cell
  group_by(Species, Haplotype, Chr) %>%
  slice_max(order_by = Length, n = 1, with_ties = FALSE) %>%
  ungroup()

# 4. Define Strict Factor Ordering Logic
# Define Species order
species_levels <- c(
  "O. sativa ssp. japonica", "O. rufipogon", "O. sativa ssp. indica", 
  "O. nivara", "O. glaberrima", "O. glumaepatula", "O. longistaminata", 
  "O. punctata", "O. officinalis", "O. australiensis", "O. brachyantha", 
  "O. meyeriana", "L. hexandra"
)

# Define Chromosome order (Chr01-Chr12)
chr_levels <- paste0("Chr", sprintf("%02d", 1:12))

# Define Haplotype order: Levels set from 1 to 4; 
# integrated with facet_grid to ensure hap1 is positioned at the top.
hap_levels <- c("hap1", "hap2", "hap3", "hap4", "-")

df_final <- df_clean %>%
  mutate(
    Species = factor(Species, levels = species_levels),
    Chr = factor(Chr, levels = chr_levels),
    Haplotype = factor(Haplotype, levels = hap_levels)
  )

# 5. Data Visualization
p <- ggplot(df_final, aes(x = Chr, y = 1, fill = Length)) +
  geom_tile(color = "white", size = 0.5) +
  # Use Species + Haplotype faceting to achieve physical row isolation
  # switch = "y" places labels on the left side of the plot
  facet_grid(Species + Haplotype ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(
    colors = c("#313695", "#4575b4", "#abd9e9", "#fee090", "#f46d43", "#a50026"),
    name = "Length (bp)"
  ) +
  theme_minimal() +
  labs(x = "Chromosome", y = NULL) +
  theme(
    # Adjust left-side species and haplotype text: non-rotated, italicized, and right-aligned
    strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1, size = 9),
    # Hide primary y-axis labels and ticks
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # Rotate chromosome labels by 45 degrees for better readability
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Crucial: Control the vertical spacing between haplotype panels
    panel.spacing.y = unit(0.3, "lines"), 
    panel.grid = element_blank(),
    # Facet strip background and placement configuration
    strip.background = element_rect(fill = "gray98", color = "white"),
    strip.placement = "outside"
  )

# 6. Execute and Preview Plot
print(p)

# 7. Export Visualization (Large height is recommended due to the high number of facets)
ggsave("Rice_Centromere_Heatmap_Final.pdf", p, width = 11, height = 20)

--- END OF FILE ---