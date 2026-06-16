#!/usr/bin/env Rscript

# Load necessary libraries
library(optparse)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Set command-line arguments
option_list <- list(
  make_option(c("--element"), type="character", default=NULL, help="Element composition file (TSV)", metavar="character"),
  make_option(c("--cen"), type="character", default=NULL, help="Centromere region coordinates file (TSV)", metavar="character"),
  make_option(c("--output"), type="character", default=NULL, help="Output figure filename (e.g., output.pdf or output.png)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$element) || is.null(opt$cen) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("All input and output parameters must be provided.", call.=FALSE)
}

# 2. Define color palette
color_dict <- c(
  'Satellite'     = '#D62728',
  'Intact LTR-RT' = '#1F77B4',
  'NUMT'          = '#FFBB78',
  'NUPT'          = '#8C564B',
  '5S rDNA'       = '#2CA02C',
  'Gene'          = '#9467BD',
  'Other'         = '#C7C7C7'
)

# Define stacking order for the bar plot
stack_order <- c('Satellite', 'Intact LTR-RT', 'NUMT', 'NUPT', '5S rDNA', 'Gene', 'Other')

# Define species order (fixed sequence)
species_order <- c(
  "O. sativa ssp. japonica", "O. rufipogon_hap1", "O. rufipogon_hap2",
  "O. sativa ssp. indica", "O. nivara_hap1", "O. nivara_hap2",
  "O. glaberrima_hap1", "O. glaberrima_hap2", "O. glumaepatula_hap1",
  "O. glumaepatula_hap2", "O. longistaminata_hap1", "O. longistaminata_hap2",
  "O. punctata_hap1", "O. punctata_hap2", "O. officinalis_hap1",
  "O. officinalis_hap2", "O. australiensis_hap1", "O. australiensis_hap2",
  "O. brachyantha_hap1", "O. brachyantha_hap2", "O. meyeriana_hap1",
  "O. meyeriana_hap2", "L. hexandra_hap1", "L. hexandra_hap3",
  "L. hexandra_hap2", "L. hexandra_hap4"
)

# 3. Read and process data
df_cen <- read.table(opt$cen, sep="\t", header=FALSE, 
                     col.names=c("Chr", "Start", "End", "assembly"), stringsAsFactors=FALSE)

df_total <- df_cen %>%
  mutate(Region_Len = End - Start) %>%
  group_by(assembly) %>%
  summarise(Total_Len = sum(Region_Len))

df_el <- read.table(opt$element, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

# Function for element type mapping
map_type <- function(t) {
  t_low <- tolower(as.character(t))
  if (t_low == "satellite") return("Satellite")
  if (t_low == "intactltr") return("Intact LTR-RT")
  if (t_low == "numt") return("NUMT")
  if (t_low == "nupt") return("NUPT")
  if (grepl("5s", t_low)) return("5S rDNA")
  if (t_low == "gene") return("Gene")
  return(NA)
}

# Summarize identified elements
df_el_clean <- df_el %>%
  mutate(Type = sapply(ElementType, map_type)) %>%
  filter(!is.na(Type)) %>%
  group_by(assembly, Type) %>%
  summarise(Length = sum(LengthInCen), .groups = 'drop')

# 4. Calculate "Other" sequences (unidentified or remaining portions)
df_identified_sum <- df_el_clean %>%
  group_by(assembly) %>%
  summarise(Identified_Sum = sum(Length))

df_other <- df_total %>%
  left_join(df_identified_sum, by="assembly") %>%
  replace_na(list(Identified_Sum = 0)) %>%
  mutate(Length = Total_Len - Identified_Sum) %>%
  mutate(Length = ifelse(Length < 0, 0, Length)) %>%
  mutate(Type = "Other") %>%
  select(assembly, Type, Length)

df_plot <- bind_rows(df_el_clean, df_other) %>%
  mutate(Length_Mb = Length / 1e6)

# 5. Factor levels and order control
df_plot$assembly <- factor(df_plot$assembly, levels = rev(species_order))
df_plot$Type <- factor(df_plot$Type, levels = stack_order)

# 6. Data visualization
p <- ggplot(df_plot, aes(x = assembly, y = Length_Mb, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, linewidth = 0.1) +
  scale_fill_manual(values = color_dict, name = "Sequence type") +
  coord_flip() +
  theme_classic() +
  labs(x = "", y = "Centromere length (Mb)") +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11, margin = margin(t = 10)),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  )

# 7. Save output file
num_species <- length(unique(df_plot$assembly))
fig_height <- max(4, num_species * 0.3)

ggsave(opt$output, plot = p, width = 7, height = fig_height, dpi = 300)

cat(paste0("Plotting successful! Result saved to: ", opt$output, "\n"))