file_path <- "/share/org/YZWL/yzwl_shahd/YanRY/CENH3_bdg_10b_raw/CC_EE_sat/CEN155_EE_CC_CEN126_CC.csv"

data <- read.csv(
  file_path,
  header = TRUE,
  colClasses = c("numeric", "numeric", "numeric")
)

# Verify data loading integrity
head(data)
summary(data)
library(dplyr)
library(ggplot2)

# Reshape data to long format
data_long <- data %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  mutate(
    group = ifelse(Sample == "CEN155_EE", "EE_CEN155",
                   ifelse(Sample == "CEN155_CC", "CC_CEN155", "CC_CEN126"))
  )

# Generate violin plots (with embedded boxplots)
p <- ggplot(data_long, aes(x = group, y = Value, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  labs(
    x = "Group",
    y = "Average Signal",
    title = "Average Signal Distribution by Group",
    fill = "group"
  ) +
  scale_fill_manual(values = c(
    "CC_CEN126" = "#84ba42",  # Green
    "CC_CEN155" = "#b85c5c",  # Red
    "EE_CEN155" = "#d4a06f"   # Brown
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "right",
    legend.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# Display plot
print(p)