library(ggplot2)
library(dplyr)

# Data loading and preparation
cenh3 <- read.csv("fai_summary.csv", header = TRUE) %>%
  mutate(
    len_mb = L/S,  
    Chromosome = factor(Chromosome, levels = unique(Chromosome)),
    Species_group = gsub("_hap[12]$", "", Species)  # Remove the _hap1/_hap2 suffixes
  )

# Define species groups and color mapping
species_groups <- c("AA_Oruf", "AA_Ogla", "AA_Oniv", "AA_Olon", 
                    "AA_Oglu", "BB_Opun", "CC_Ooff", "EE_Oaus",
                    "FF_Obra", "GG_Omey", "XX_Lhex")

color_mapping <- c(
  'AA_Oruf' = "#215A20",
  'AA_Ogla' = "#76D273",
  'AA_Oniv' = "#3BA738",
  'AA_Olon' = "#51C54E",
  'AA_Oglu' = "#3D8347",
  'BB_Opun' = "#F2AE2C",
  'CC_Ooff' = "#684E94",
  'EE_Oaus' = "#4E84C3",
  'FF_Obra' = "#D55F6F",
  'GG_Omey' = "#9D5427",
  'XX_Lhex' = "#595959"
)

# Plotting
ggplot(cenh3, aes(x = Chromosome, y = len_mb)) +
  # Generate violin plots
  geom_violin(scale = "width", trim = TRUE, color = "gray6", alpha = 0.7) +
  
  # Add scatter points stratified by species
  geom_point(aes(color = Species_group), 
             position = position_jitter(width = 0.2, height = 0), 
             size = 2, shape = 15) +
  
  # Apply the color scheme
  scale_color_manual(values = color_mapping,
                     breaks = names(color_mapping),
                     labels = names(color_mapping)) +
  
  # Labels and themes
  labs(x = "Chromosome", y = "L/S", color = "Species Group") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.5, "cm")
  ) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))