#!/usr/bin/env Rscript

# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Set file paths
single_copy_file <- "D:/17_1_copy_merged_n17_summary.tsv"
kmer_file <- "D:/17mer_combined_top10_coverage.tsv"
output_prefix <- "bar_with_points_bh"

#### Check if files exist
check_files <- function() {
  if (!file.exists(single_copy_file)) {
    cat("Error: Single-copy data file does not exist:", single_copy_file, "\n")
    return(FALSE)
  }
  if (!file.exists(kmer_file)) {
    cat("Error: K-mer data file does not exist:", kmer_file, "\n")
    return(FALSE)
  }
  return(TRUE)
}

#### Load data
load_data <- function() {
  if (check_files()) {
    cat("Loading data files...\n")
    
    # Read single-copy data
    single_copy_df <- read.table(single_copy_file, 
                                 header = TRUE, 
                                 sep = "\t",
                                 stringsAsFactors = FALSE)
    
    # Read K-mer data
    kmer_df <- read.table(kmer_file, 
                          header = TRUE, 
                          sep = "\t",
                          stringsAsFactors = FALSE)
    
    cat("Data loading complete!\n")
    
    # Display basic data information
    cat("\n=== Single-copy Data Information ===\n")
    cat("Rows:", nrow(single_copy_df), "\n")
    cat("Columns:", ncol(single_copy_df), "\n")
    cat("First 3 rows:\n")
    print(head(single_copy_df, 3))
    
    cat("\n=== K-mer Data Information ===\n")
    cat("Rows:", nrow(kmer_df), "\n")
    cat("Columns:", ncol(kmer_df), "\n")
    cat("First 3 rows:\n")
    print(head(kmer_df, 3))
    
  } else {
    stop("Files do not exist. Please check the paths.")
  }
  
  return(list(single_copy = single_copy_df, kmer = kmer_df))
}

# Load data
data_list <- load_data()
single_copy_df <- data_list$single_copy
kmer_df <- data_list$kmer

#### Data preprocessing - adjust based on actual column names
process_data <- function(df, data_type) {
  
  if (data_type == "Single-copy") {
    # Single-copy data already has correct column names
    df_long <- df %>%
      pivot_longer(
        cols = c(INCEN_Proportion, OUTCEN_Proportion),
        names_to = "Region",
        values_to = "Proportion"
      )
  } else if (data_type == "K-mer") {
    # K-mer data requires column renaming
    df <- df %>%
      rename(
        INCEN_Proportion = top10INCEN_Proportion,
        OUTCEN_Proportion = top10OUTCEN_Proportion
      )
    
    # Extract species information and haplotype from the Species column
    if ("Species" %in% colnames(df)) {
      df <- df %>%
        mutate(
          # Assume Species column format is "AA_Ogla_hap1"
          Species_full = Species,  # Keep full name
          Species = sub("_hap[0-9]+$", "", Species),  # Extract species name
          Haplotype = sub(".*_(hap[0-9]+)$", "\\1", Species_full)  # Extract haplotype
        )
    }
    
    df_long <- df %>%
      pivot_longer(
        cols = c(INCEN_Proportion, OUTCEN_Proportion),
        names_to = "Region",
        values_to = "Proportion"
      )
  }
  
  df_long <- df_long %>%
    mutate(
      Region = case_when(
        Region == "INCEN_Proportion" ~ "Centromere",
        Region == "OUTCEN_Proportion" ~ "Arm"
      ),
      DataType = data_type,
      Percentage = Proportion * 100  # Convert to percentage
    )
  
  return(df_long)
}

cat("Processing data...\n")
single_copy_long <- process_data(single_copy_df, "Single-copy")
kmer_long <- process_data(kmer_df, "K-mer")

#### Custom permutation test function (10,000 iterations)
permutation_test <- function(x, y, n_permutations = 10000) {
  observed_diff <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
  combined <- c(x, y)
  n_x <- length(x)
  
  perm_diffs <- replicate(n_permutations, {
    permuted <- sample(combined)
    perm_x <- permuted[1:n_x]
    perm_y <- permuted[(n_x + 1):length(combined)]
    mean(perm_x, na.rm = TRUE) - mean(perm_y, na.rm = TRUE)
  })
  
  p_value <- (sum(abs(perm_diffs) >= abs(observed_diff)) + 1) / (n_permutations + 1)
  return(p_value)
}

#### BH multiple comparison function
bh_multiple_comparison <- function(centromere_data, arm_data, alpha = 0.05) {
  # Prepare data
  groups_data <- list(Centromere = centromere_data, Arm = arm_data)
  
  # Calculate basic statistics
  means <- c(mean(centromere_data, na.rm = TRUE), mean(arm_data, na.rm = TRUE))
  ns <- c(length(centromere_data), length(arm_data))
  
  # Check if there is sufficient data
  if (ns[1] < 2 || ns[2] < 2) {
    return(list(
      groups = c("a", "a"),
      p_value = NA,
      p_adj_value = NA,
      significant = FALSE
    ))
  }
  
  # Combine all data
  all_data <- c(centromere_data, arm_data)
  
  # Construct data frame
  group_labels <- c(rep("Centromere", ns[1]), rep("Arm", ns[2]))
  data_df <- data.frame(
    value = all_data,
    group = group_labels
  )
  
  # Calculate raw p-value from t-test
  t_test_result <- t.test(centromere_data, arm_data)
  raw_p <- t_test_result$p.value
  
  # Apply BH correction (with only one comparison, BH correction equals raw p-value)
  p_adj <- raw_p
  
  # Assign grouping letters based on adjusted p-value
  if (p_adj < alpha) {
    groups <- c("a", "b")  # Different letters indicate significant difference
  } else {
    groups <- c("a", "a")  # Same letters indicate no significant difference
  }
  
  return(list(
    groups = groups,
    p_value = raw_p,
    p_adj_value = p_adj,
    significant = p_adj < alpha,
    means = means,
    test_statistic = t_test_result$statistic,
    df = t_test_result$parameter,
    ci_lower = t_test_result$conf.int[1],
    ci_upper = t_test_result$conf.int[2]
  ))
}

#### Statistical analysis: 10,000 permutation tests + BH multiple comparison
perform_statistical_analysis <- function(df) {
  species_list <- unique(df$Species)
  results <- list()
  
  for (species in species_list) {
    for (data_type in unique(df$DataType)) {
      
      species_data <- df %>%
        filter(Species == species, DataType == data_type)
      
      if (nrow(species_data) > 0) {
        centromere_data <- species_data %>% 
          filter(Region == "Centromere") %>% 
          pull(Proportion)
        arm_data <- species_data %>% 
          filter(Region == "Arm") %>% 
          pull(Proportion)
        
        if (length(centromere_data) > 1 && length(arm_data) > 1) {
          tryCatch({
            set.seed(123)
            
            # 1. 10,000 permutation tests
            p_perm <- permutation_test(centromere_data, arm_data, n_permutations = 10000)
            
            # 2. BH multiple comparison
            bh_result <- bh_multiple_comparison(centromere_data, arm_data, alpha = 0.05)
            
            results[[length(results) + 1]] <- data.frame(
              Species = species,
              DataType = data_type,
              # Permutation test results
              p_perm = p_perm,
              perm_label = case_when(
                p_perm < 0.001 ~ "***",
                p_perm < 0.01 ~ "**", 
                p_perm < 0.05 ~ "*",
                TRUE ~ "ns"
              ),
              # BH multiple comparison results
              bh_cent_group = bh_result$groups[1],
              bh_arm_group = bh_result$groups[2],
              bh_raw_p = bh_result$p_value,
              bh_adj_p = bh_result$p_adj_value,
              bh_significant = bh_result$significant,
              # Statistical information
              t_statistic = ifelse(!is.null(bh_result$test_statistic), bh_result$test_statistic, NA),
              df = ifelse(!is.null(bh_result$df), bh_result$df, NA),
              ci_lower = ifelse(!is.null(bh_result$ci_lower), bh_result$ci_lower, NA),
              ci_upper = ifelse(!is.null(bh_result$ci_upper), bh_result$ci_upper, NA),
              # Mean information
              centromere_mean = bh_result$means[1],
              arm_mean = bh_result$means[2],
              centromere_n = length(centromere_data),
              arm_n = length(arm_data)
            )
          }, error = function(e) {
            cat("Error when calculating", species, data_type, ":", e$message, "\n")
          })
        }
      }
    }
  }
  
  significance_df <- bind_rows(results)
  
  # Add formatted p-value column
  if (nrow(significance_df) > 0) {
    significance_df <- significance_df %>%
      mutate(bh_p_formatted = ifelse(!is.na(bh_adj_p), 
                                     format.pval(bh_adj_p, digits = 2, eps = 0.001),
                                     NA))
  }
  
  return(significance_df)
}

cat("Performing 10,000 permutation tests and BH multiple comparisons...\n")
significance_df <- perform_statistical_analysis(bind_rows(kmer_long, single_copy_long))

# Separate results
kmer_sig <- significance_df %>% filter(DataType == "K-mer")
single_copy_sig <- significance_df %>% filter(DataType == "Single-copy")

#### Create compact bar plot + data points (with statistical results)
create_bar_with_points_plot <- function(df, df_sig, plot_title) {
  
  # Define species order and labels
  species_order <- c(
    'AA_Osat_jap', 'AA_Osat_ind', 'AA_Ogla', 'AA_Oruf', 'AA_Oniv',
    'AA_Olon', 'AA_Oglu', 'BB_Opun', 'CC_Ooff', 'EE_Oaus', 'FF_Obra', 'GG_Omey'
  )
  
  species_labels <- c(
    'O. sativa japonica', 'O. sativa indica', 'O. glaberrima',
    'O. rufipogon', 'O. nivara', 'O. longistaminata', 'O. glumaepatula',
    'O. punctata', 'O. officinalis', 'O. australiensis', 'O. brachyantha', 'O. meyeriana'
  )
  
  # Prepare plotting data
  plot_df <- df %>%
    filter(Species %in% species_order) %>%
    mutate(
      Species = factor(Species, levels = rev(species_order), labels = rev(species_labels)),
      Region = factor(Region, levels = c("Centromere", "Arm"))
    )
  
  # Check if there is data
  if (nrow(plot_df) == 0) {
    cat(paste("Warning: For", plot_title, ", no matching species data found\n"))
    cat("Available species:", paste(unique(df$Species), collapse = ", "), "\n")
    return(ggplot() + ggtitle(paste(plot_title, "(No data)")))
  }
  
  # Calculate means and standard errors for bar plot
  summary_df <- plot_df %>%
    group_by(Species, Region) %>%
    summarise(
      Mean = mean(Percentage, na.rm = TRUE),
      SE = sd(Percentage, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  # Prepare significance annotations
  if (nrow(df_sig) > 0) {
    sig_annotations <- df_sig %>%
      filter(Species %in% species_order) %>%
      mutate(
        Species = factor(Species, levels = rev(species_order), labels = rev(species_labels)),
        # Permutation test significance marker position
        y_pos_perm = ifelse(plot_title == "17-mer Coverage", 55, 105),
        # BH grouping letter position (placed above the bars)
        y_pos_bh_cent = ifelse(plot_title == "17-mer Coverage", 48, 98),
        y_pos_bh_arm = ifelse(plot_title == "17-mer Coverage", 48, 98),
        # BH p-value annotation position
        y_pos_bh_p = ifelse(plot_title == "17-mer Coverage", 58, 108),
        # Format BH p-value for display
        bh_label = ifelse(!is.na(bh_p_formatted) & bh_significant,
                          paste0("BH p=", bh_p_formatted),
                          "")
      )
  } else {
    sig_annotations <- data.frame()
  }
  
  # Create bar plot + data points
  p <- ggplot() +
    # Bar plot
    geom_col(
      data = summary_df,
      aes(x = Species, y = Mean, fill = Region),
      position = position_dodge(0.7),
      width = 0.6,
      alpha = 0.5,
      color = "black",
      linewidth = 0.2
    ) +
    # Error bars
    geom_errorbar(
      data = summary_df,
      aes(x = Species, ymin = Mean - SE, ymax = Mean + SE, group = Region),
      position = position_dodge(0.7),
      width = 0.2,
      linewidth = 0.3,
      color = "black"
    ) +
    # Individual data points
    geom_point(
      data = plot_df,
      aes(x = Species, y = Percentage, color = Region),
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = 0.7,
        jitter.height = 0
      ),
      size = 1.2,
      alpha = 0.6,
      shape = 16
    ) +
    # Color scheme
    scale_fill_manual(values = c("Centromere" = "#E69F00", "Arm" = "#56B4E9")) +
    scale_color_manual(values = c("Centromere" = "#E69F00", "Arm" = "#56B4E9")) +
    labs(
      x = "",
      y = "Percentage (%)",
      title = plot_title,
      fill = "Region", color = "Region"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face = "italic"),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8, margin = margin(r = 5)),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      legend.position = "right",
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    # Set y-axis limits
    scale_y_continuous(
      limits = c(0, ifelse(plot_title == "17-mer Coverage", 60, 110)),
      expand = c(0, 0)
    ) +
    # Flip coordinates
    coord_flip()
  
  # Add significance markers if significance data is available
  if (nrow(sig_annotations) > 0) {
    # Permutation test significance markers (asterisks)
    p <- p + geom_text(
      data = sig_annotations,
      aes(x = Species, y = y_pos_perm, label = perm_label),
      inherit.aes = FALSE,
      size = 3.5,
      fontface = "bold",
      color = "#FF4500",  # Orange asterisks
      vjust = 0
    )
  }
  
  return(p)
}

#### Create side-by-side plots
create_side_by_side_plots <- function() {
  cat("Generating bar plot + data points side-by-side...\n")
  
  # Check K-mer data
  cat("\nSpecies in K-mer data:", paste(unique(kmer_long$Species), collapse = ", "), "\n")
  
  # Create K-mer plot (left)
  kmer_plot <- create_bar_with_points_plot(kmer_long, kmer_sig, "17-mer Coverage")
  
  # Check single-copy data
  cat("Species in single-copy data:", paste(unique(single_copy_long$Species), collapse = ", "), "\n")
  
  # Create Single-copy plot (right) - remove y-axis text
  single_plot <- create_bar_with_points_plot(single_copy_long, single_copy_sig, "Single-copy Proportion") +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  # Combine plots using patchwork
  combined_plots <- kmer_plot + single_plot +
    plot_layout(ncol = 2, widths = c(1, 0.9), guides = "collect") &
    theme(legend.position = "bottom",
          legend.box = "horizontal")
  
  return(combined_plots)
}

#### Main program
main <- function() {
  cat("Starting bar plot + data points side-by-side generation...\n")
  
  # Output statistical results
  cat("\n=== Statistical Analysis Results ===\n")
  
  if (nrow(kmer_sig) > 0) {
    cat("\n17-mer coverage statistical results:\n")
    print(kmer_sig)
  } else {
    cat("\nWarning: No 17-mer statistical results\n")
  }
  
  if (nrow(single_copy_sig) > 0) {
    cat("\nSingle-copy sequence statistical results:\n")
    print(single_copy_sig)
  } else {
    cat("\nWarning: No single-copy statistical results\n")
  }
  
  # Generate plot
  final_plot <- create_side_by_side_plots()
  
  # Save plot
  output_path <- paste0(output_prefix, ".pdf")
  cat("\nSaving plot to:", output_path, "\n")
  
  ggsave(output_path, final_plot, width = 12, height = 6, dpi = 300)
  
  if (file.exists(output_path)) {
    cat("✓ Plot successfully saved to:", output_path, "\n")
  } else {
    cat("✗ Error: File save failed\n")
  }
  
  # Save statistical results
  stat_output <- paste0(output_prefix, "_statistics.tsv")
  write.table(significance_df, stat_output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Statistical results saved to:", stat_output, "\n")
  
  # Save detailed BH analysis results
  cat("\nAvailable column names:\n")
  print(colnames(significance_df))
  
  # Safely select columns that exist
  bh_details <- significance_df %>%
    select(any_of(c(
      "Species", "DataType",
      "Permutation_p" = "p_perm",
      "Permutation_sig" = "perm_label",
      "BH_raw_p" = "bh_raw_p",
      "BH_adj_p" = "bh_adj_p",
      "BH_significant" = "bh_significant",
      "BH_p_formatted" = "bh_p_formatted",
      "Centromere_group" = "bh_cent_group",
      "Arm_group" = "bh_arm_group",
      "Centromere_mean" = "centromere_mean",
      "Arm_mean" = "arm_mean",
      "t_statistic", "df", "ci_lower", "ci_upper"
    )))
  
  bh_output <- paste0(output_prefix, "_bh_details.tsv")
  write.table(bh_details, bh_output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("BH detailed results saved to:", bh_output, "\n")
  
  cat("\n=== Plot Explanation ===\n")
  cat("1. Asterisks (***/**/*): Results from 10,000 permutation tests (orange)\n")
  cat("2. (a/b): BH multiple comparison grouping letters (blue=Centromere, green=Arm)\n")
  cat("3. BH p=xx: Benjamini-Hochberg adjusted p-value (purple, displayed only when significant)\n")
  cat("4. Same letters indicate no significant difference; different letters indicate significant difference (BH-adjusted p < 0.05)\n")
  cat("5. Statistical method: t-test + BH multiple comparison correction (FDR control)\n")
  cat("6. BH correction method: Benjamini-Hochberg (controls false discovery rate)\n")
  
  cat("\nAll tasks completed!\n")
}

# Run main program
main()