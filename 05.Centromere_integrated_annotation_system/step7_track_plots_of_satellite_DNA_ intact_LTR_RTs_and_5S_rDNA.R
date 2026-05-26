# ==============================================================================
# Genomic Track Plot Generator for BB_Opun (V2 - Centromere Highlight)
#
# Version Notes:
# - New: Added semi-transparent grey rectangles to the background of each track 
#   to highlight centromere core regions.
# - The script automatically retrieves centromere coordinates for the current 
#   chromosome from the predefined list.
# ==============================================================================

# ==============================================================================
# 1. Load Required Libraries
# ==============================================================================
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(stringr)

# ==============================================================================
# 2. Parameter Configuration
# ==============================================================================
# !!! [Please modify according to your specific environment] !!!
# Set the local root directory containing all input files
base_path <- "C:/Users/DELL/Desktop/轨道/BB"   # <--- Assumes files and script are in the same directory

# Set species name prefix for file identification
species_prefix <- "BB_Opun"

# PDF output parameters
pdf_width <- 30
pdf_height <- 12

# --- Define all regions to be plotted ---
regions_to_plot <- data.frame(
  hap = c("hap1", "hap1"),
  chr = c("Chr11", "Chr10"),
  start = c(14012000, 10440000),
  end = c(16442000, 11746000)
)

# --- New: Define centromere core regions for all chromosomes ---
cen_regions_to_highlight <- data.frame(
  hap = c("hap1", "hap2", "hap1", "hap2", "hap1", "hap2", "hap1", "hap2", "hap1", "hap2", "hap1", "hap2", "hap1", "hap2", "hap1", "hap2", "hap1", "hap2", "hap1", "hap2", "hap1", "hap2", "hap1", "hap2"),
  chr = rep(paste0("Chr", sprintf("%02d", 1:12)), each = 2),
  start = c(18498000, 18512000, 16200000, 16206000, 22442000, 22414000, 11508000, 11686000, 13010000, 13016000, 20010000, 20006000, 13490000, 13490000, 13728000, 13726000, 5430000, 5458000, 10940000, 10808000, 14512000, 14548000, 12928000, 12918000),
  end = c(19218000, 19194000, 16960000, 16974000, 23252000, 23264000, 12274000, 12466000, 13628000, 13628000, 20498000, 20498000, 13920000, 13922000, 14292000, 14304000, 6038000, 6062000, 11246000, 11110000, 15942000, 15912000, 13618000, 13656000)
)


# ==============================================================================
# 3. Data Loading (Consistent with previous versions)
# ==============================================================================
cat("Loading data for", species_prefix, "...\n")
cenh3_files <- list.files(base_path, pattern = ".CENH3.bdg$", full.names = TRUE); all_cenH3 <- lapply(cenh3_files, function(f) { hap_type <- str_extract(basename(f), "hap[12]"); rep_type <- paste0("rep", str_extract(basename(f), "(?<=sample)[123]")); read.table(f, header = FALSE, col.names = c("chr", "start", "end", "log2ratio")) %>% mutate(hap = hap_type, replicate = rep_type) }) %>% bind_rows(); cat("CENH3 data loading completed.\n")
sat_files <- list.files(base_path, pattern = paste0(species_prefix, ".*.sat.INcen.bed"), full.names = TRUE); all_satellites <- lapply(sat_files, function(f) { hap_type <- ifelse(grepl("hap1", f), "hap1", "hap2"); read.table(f, header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "id", "name", "strand")) %>% mutate(hap = hap_type) }) %>% bind_rows(); cat("Satellite data loading completed.\n")
rdna_files <- list.files(base_path, pattern = paste0(species_prefix, ".*.rDNA.bed"), full.names = TRUE); all_rdna <- lapply(rdna_files, function(f) { hap_type <- ifelse(grepl("hap1", f), "hap1", "hap2"); read.table(f, header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand")) %>% mutate(hap = hap_type) }) %>% bind_rows(); cat("5S rDNA data loading completed.\n")
ltr_files <- list.files(base_path, pattern = paste0(species_prefix, ".*.intactLTR.bed"), full.names = TRUE); all_ltrs <- lapply(ltr_files, function(f) { hap_type <- ifelse(grepl("hap1", f), "hap1", "hap2"); read.table(f, header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "te_type", "element_id", "strand", "attributes"), quote = "", comment.char = "") %>% mutate(hap = hap_type) %>% mutate(tsd = str_match(attributes, "tsd=([A-Z]+)")[, 2]) }) %>% bind_rows(); all_ltrs$tsd[is.na(all_ltrs$tsd)] <- "unknown"; cat("Intact TEs data loading completed.\n")
hite_files <- list.files(base_path, pattern = paste0(species_prefix, ".*.HiTE.gff"), full.names = TRUE); all_hite <- lapply(hite_files, function(f) { hap_type <- ifelse(grepl("hap1", f), "hap1", "hap2"); read.table(f, header = FALSE, sep = "\t", comment.char = "#", col.names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")) %>% mutate(hap = hap_type) %>% select(chr, start, end, hap) }) %>% bind_rows(); cat("Fragment TEs data loading completed.\n")
cat("All data loading tasks completed.\n")


# ==============================================================================
# 4. Plotting Loop
# ==============================================================================
for (i in 1:nrow(regions_to_plot)) {
  
  current_hap <- regions_to_plot$hap[i]
  current_chr <- regions_to_plot$chr[i]
  view_start_abs <- regions_to_plot$start[i]
  view_end_abs <- regions_to_plot$end[i]
  
  cat(paste("\nGenerating plots for region:", current_hap, current_chr, "...\n"))
  
  coord_offset <- view_start_abs
  relative_length <- view_end_abs - view_start_abs
  
  # --- New: Filter centromere regions for highlight and convert to relative coordinates ---
  cen_highlight_df <- cen_regions_to_highlight %>%
    filter(hap == current_hap, chr == current_chr) %>%
    mutate(
      rel_start = start - coord_offset,
      rel_end = end - coord_offset
    )
  
  # ... (Data filtering section remains unchanged) ...
  cenH3_filtered <- all_cenH3 %>% filter(hap == current_hap, chr == current_chr, start < view_end_abs, end > view_start_abs) %>% mutate(rel_start = start - coord_offset)
  sat_filtered <- all_satellites %>% filter(hap == current_hap, chr == current_chr, start < view_end_abs, end > view_start_abs) %>% mutate(rel_start = start - coord_offset, rel_end = end - coord_offset)
  rdna_filtered <- all_rdna %>% filter(hap == current_hap, chr == current_chr, start < view_end_abs, end > view_start_abs) %>% mutate(rel_start = start - coord_offset, rel_end = end - coord_offset)
  ltr_filtered <- all_ltrs %>% filter(hap == current_hap, chr == current_chr, start < view_end_abs, end > view_start_abs) %>% mutate(rel_start = start - coord_offset, rel_end = end - coord_offset)
  hite_filtered <- all_hite %>% filter(hap == current_hap, chr == current_chr, start < view_end_abs, end > view_start_abs) %>% mutate(rel_start = start - coord_offset, rel_end = end - coord_offset)
  
  # ... (Color and theme definitions remain unchanged) ...
  theme_track <- theme_minimal() + theme(panel.border = element_rect(colour = "black", fill=NA), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_blank(), legend.position = "right", axis.title.y = element_text(size = 14, angle = 0, vjust = 0.5, hjust=1))
  unique_ltr_types <- unique(all_ltrs$te_type); ltr_colors <- colorRampPalette(brewer.pal(min(8, length(unique_ltr_types)), "Set2"))(length(unique_ltr_types)); names(ltr_colors) <- unique_ltr_types
  unique_tsd_types <- unique(all_ltrs$tsd); tsd_colors <- colorRampPalette(brewer.pal(min(12, length(unique_tsd_types)), "Paired"))(length(unique_tsd_types)); names(tsd_colors) <- unique_tsd_types
  
  # --- Revised: Add background highlight layer to each track plot ---
  p_cenh3 <- ggplot() + geom_rect(data = cen_highlight_df, aes(xmin = rel_start, xmax = rel_end, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.5) + geom_line(data=cenH3_filtered, aes(x = rel_start, y = log2ratio, color = replicate), size = 0.6) + scale_color_manual(values = c("rep1" = "#008C8C", "rep2" = "#dbb428", "rep3" = "#d4562e"), name = "CENH3 Replicate") + scale_x_continuous(limits = c(0, relative_length), expand = c(0, 0)) + ylab("CENH3") + theme_track
  p_sat <- ggplot() + geom_rect(data = cen_highlight_df, aes(xmin = rel_start, xmax = rel_end, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.5) + geom_rect(data=sat_filtered, aes(xmin = rel_start, xmax = rel_end, ymin = -0.5, ymax = 0.5, fill = name)) + scale_x_continuous(limits = c(0, relative_length), expand = c(0, 0)) + scale_y_continuous(limits = c(-1, 1), breaks = NULL) + scale_fill_brewer(palette = "Set1", name = "Satellite Type") + ylab("Satellite") + theme_track
  p_rdna <- ggplot() + geom_rect(data = cen_highlight_df, aes(xmin = rel_start, xmax = rel_end, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.5) + geom_rect(data=rdna_filtered, aes(xmin = rel_start, xmax = rel_end, ymin = -0.5, ymax = 0.5, fill = name)) + scale_x_continuous(limits = c(0, relative_length), expand = c(0, 0)) + scale_y_continuous(limits = c(-1, 1), breaks = NULL) + scale_fill_brewer(palette = "Accent", name = "rDNA Type") + ylab("5S rDNA") + theme_track
  p_ltr_intact <- ggplot() + geom_rect(data = cen_highlight_df, aes(xmin = rel_start, xmax = rel_end, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.5) + geom_rect(data=ltr_filtered, aes(xmin = rel_start, xmax = rel_end, ymin = -0.5, ymax = 0.5, fill = te_type)) + scale_x_continuous(limits = c(0, relative_length), expand = c(0, 0)) + scale_y_continuous(limits = c(-1, 1), breaks = NULL) + scale_fill_manual(values = ltr_colors, name = "Intact TE Superfamily") + ylab("Intact TEs") + theme_track
  p_tsd <- ggplot() + geom_rect(data = cen_highlight_df, aes(xmin = rel_start, xmax = rel_end, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.5) + geom_rect(data=ltr_filtered, aes(xmin = rel_start, xmax = rel_end, ymin = -0.5, ymax = 0.5, fill = tsd)) + scale_x_continuous(limits = c(0, relative_length), expand = c(0, 0)) + scale_y_continuous(limits = c(-1, 1), breaks = NULL) + scale_fill_manual(values = tsd_colors, name = "TSD Sequence") + ylab("TSD") + theme_track
  p_frag_tes <- ggplot() + geom_rect(data = cen_highlight_df, aes(xmin = rel_start, xmax = rel_end, ymin = -Inf, ymax = Inf), fill = "lightgrey", alpha = 0.5) + geom_rect(data=hite_filtered, aes(xmin = rel_start, xmax = rel_end, ymin = -0.5, ymax = 0.5), fill = "grey50") + scale_x_continuous(limits = c(0, relative_length), expand = c(0, 0)) + scale_y_continuous(limits = c(-1, 1), breaks = NULL) + ylab("Fragment\nTEs") + theme_track
  
  # ... (X-axis addition, combination, and saving sections remain unchanged) ...
  x_axis_breaks <- seq(from = 0, to = relative_length, by = 200 * 1000)
  x_axis_labels <- (view_start_abs + x_axis_breaks) / 1000000
  p_frag_tes_with_axis <- p_frag_tes +
    scale_x_continuous(name = paste("Position on", current_chr, "(Mb)"), limits = c(0, relative_length), expand = c(0, 0), breaks = x_axis_breaks, labels = sprintf("%.2f", x_axis_labels)) +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.ticks.x = element_line(), axis.title.x = element_text(size = 14))
  combined_plot <- (p_cenh3 / p_sat / p_rdna / p_ltr_intact / p_tsd / p_frag_tes_with_axis) +
    plot_layout(heights = c(3, 1, 1, 1, 1, 1.5), guides = 'collect') &
    theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10), legend.key.size = unit(0.5, "cm"))
  combined_plot <- combined_plot + plot_annotation(title = paste("Genomic Tracks for", species_prefix, current_hap, current_chr), theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold")))
  output_filename <- paste0(species_prefix, "_", current_hap, "_", current_chr, "_track_plot.pdf")
  ggsave(
    output_filename,
    plot = combined_plot,
    width = pdf_width,
    height = pdf_height,
    units = "in",
    device = "pdf"
  )
  cat(paste("Plot saved as:", output_filename, "\n"))
}

cat("\nAll PDF plots generated successfully!\n")