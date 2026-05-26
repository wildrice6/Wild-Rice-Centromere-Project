# Load necessary libraries
library(ggplot2)
library(dplyr) # Using dplyr for efficient data manipulation

# Step 1: Add a grouping identifier column to each dataset
outcen_identity$location <- "Outside Centromere"
incen_identity$location <- "Inside Centromere"

# Step 2: Combine the two datasets into one using rbind
combined_data <- rbind(outcen_identity, incen_identity)

# Step 3: Plotting, map color to the location column
ggplot(combined_data, aes(x = V9, y = V7)) +
    # Only a single geom_jitter layer is required now
    # Map color to location within aes()
    geom_jitter(aes(color = location), size = 0.5, alpha = 0.6, width = 0.25) +
    
    # Step 4: Use scale_color_manual() to customize colors and the legend
    scale_color_manual(
        name = "Region",  # Set the title of the legend
        values = c(
            "Outside Centromere" = "#b2bbbe", # Assign color for "Outside Centromere"
            "Inside Centromere" = "#c02c38"  # Assign color for "Inside Centromere"
        ),
        labels = c("Inside", "Outside") # Customize label text in the legend
    ) +
    
    theme_minimal() +
    coord_flip() +
    labs(x = 'Material', y = 'LTR identity')