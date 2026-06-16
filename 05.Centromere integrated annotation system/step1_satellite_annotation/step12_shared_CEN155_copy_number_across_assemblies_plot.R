# script.R
library(argparse)
library(ggplot2)
library(ggsignif)

# Create the argument parser
parser <- ArgumentParser(description="Example script")

# Add command-line arguments
parser$add_argument("--input", type="character", required=TRUE, help="Input file name")
parser$add_argument('-n', "--number", type="integer", required=TRUE, help="Number of accessions")

# Parse command-line arguments
args <- parser$parse_args()
input <- args$input
accession_number <- args$number

sharing <- read.csv(input, sep = '\t')

ggplot(sharing, aes(x=shared_across_files, y= copies_in_file))+
    geom_jitter(color='#eb4035', width=0.1, alpha=0.3, shape=21)+
    theme_bw()+
    scale_x_continuous(breaks=1:accession_number, limits=c(1,accession_number))+
    labs(x='Accessions', y='Copy numbers') -> p

output_pdf = paste0(input, '.pdf')
output_png = paste0(input, '.png')

ggsave(output_png, plot=p, width=6, height=4)
ggsave(output_pdf, plot=p, width=6, height=4)