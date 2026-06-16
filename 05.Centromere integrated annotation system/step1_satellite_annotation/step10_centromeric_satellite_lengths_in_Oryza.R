species_order <- c(
  "O. sativa" ,
  "O. glaberrima",
  "O. rufipogon",
  "O. nivara",
  "O. longistaminata",
  "O. glumaepatula",
  "O. punctata",
  "O. officinalis",
  "O. australiensis",
  "O. brachyantha",
  "O. meyeriana",
  "L. hexandra"
)

# Reverse the order using the rev() function
reversed_species_order <- rev(species_order)
# Convert the sat$V9 column into a factor with the specified level order
# The original column will be overwritten by the new factor column
sat$V9 <- factor(sat$V9, levels = reversed_species_order)

sat <- read.csv('all.sat.INcen.txt', header = FALSE, sep = '\t')

# 1. Create a named vector for color mapping
#    Converting the logic into R format
genome_colors <- c(
    'AA' = '#76D273',
    'BB' = '#F2AE2C',
    'CC' = '#684E94',
    'EE' = '#4E84C3',
    'FF' = '#D55F6F',
    'GG' = '#9D5427',
    'XX' = '#595959'
)

# 2. Modify ggplot code
ggplot(sat, aes(x = V7, y = V9, fill = V10)) + # Add fill = V10 in aes() to specify the coloring variable
    geom_density_ridges() +
    
    # Use scale_fill_manual to manually assign colors
    scale_fill_manual(
        values = genome_colors, 
        name = "Genome Type" # Optional: set legend title
    ) + 
    
    # Set labels for axes and legend
    labs(
        x = "Length (bp)",
        y = "Species"
    ) +
    
    # Use theme() to set Y-axis text to italic
    theme(
        axis.text.y = element_text(face = "italic")
    )