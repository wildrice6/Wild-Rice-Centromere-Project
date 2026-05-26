import os
import argparse
import pandas as pd

# Define command-line arguments
parser = argparse.ArgumentParser(description="Sort and process sequencing data")
parser.add_argument('--ont', type=str, required=True, help='Sequencing data file')
args = parser.parse_args()

ont = args.ont

# Input file follows the format: AA_Oruf.ont.sorted.bam

prefix = ont.split('.')[0]

# Generate BAM index
command5 = f"samtools index {prefix}.ont.sorted.bam"
os.system(command5)

# Calculate whole-genome sequencing depth/coverage
command6 = f"mosdepth -t 64 {prefix} {prefix}.ont.sorted.bam"
os.system(command6)