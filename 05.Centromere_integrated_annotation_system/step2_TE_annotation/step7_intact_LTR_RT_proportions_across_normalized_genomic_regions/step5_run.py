import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns


# Argument parsing
argparser = argparse.ArgumentParser(description='Plot LTR curves. Requires three files: 1. BED file containing LTR positional information; 2. Centromere interval file; 3. Genome size file.')
argparser.add_argument('--ltr', type=str, required=True, help='GFF3 file output by EDTA')
argparser.add_argument('--cen', type=str, required=True, help='Centromere interval file')
argparser.add_argument('--size', type=str, required=True, help='Genome size file')
args = argparser.parse_args()


def merge_size_cen(size_file, cen_file):
    """Input size and centromere files, and return a pandas DataFrame."""
    size = pd.read_csv(size_file, header=None, sep = '\t')
    size.columns = ['chr', 'size_start', 'size_end']
    cen = pd.read_csv(cen_file, header=None, sep = '\t')
    cen.columns = ['chr', 'cen_start', 'cen_end']
    size_cen = pd.merge(size, cen, on='chr')
    size_cen = size_cen.loc[:,['chr', 'size_start', 'cen_start', 'cen_end', 'size_end']]
    return size_cen

def divide_chromosome_into_bins(df: pd.DataFrame) -> pd.DataFrame:
    """
    Divide chromosomes into multiple sub-intervals (bins) according to specified rules 
    and append a sequence number to the region name.

    - size_start to cen_start: divided into 40 segments (p-arm), named L_arm_1, L_arm_2, ...
    - cen_start to cen_end: divided into 20 segments (centromere), named centromere_1, centromere_2, ...
    - cen_end to size_end: divided into 40 segments (q-arm), named R_arm_1, R_arm_2, ...

    Parameters:
    df (pd.DataFrame): DataFrame containing chromosome structural information, 
                        requiring columns 'chr', 'size_start', 'cen_start', 
                        'cen_end', and 'size_end'.

    Returns:
    pd.DataFrame: A new DataFrame where each row represents a sub-interval, 
                  containing 'chr', 'segment_start', 'segment_end', and 'region' columns.
    """
    all_bins = []  # List to store all generated sub-intervals

    # Iterate through each row of the original DataFrame
    for index, row in df.iterrows():
        chr_name = row['chr']

        # Define three regions and their corresponding number of segments
        regions = [
            # Region 1: p-arm (size_start -> cen_start)
            (row['size_start'], row['cen_start'], 40, 'L_arm'),
            # Region 2: Centromere (cen_start -> cen_end)
            (row['cen_start'], row['cen_end'], 20, 'centromere'),
            # Region 3: q-arm (cen_end -> size_end)
            (row['cen_end'], row['size_end'], 40, 'R_arm')
        ]

        # Process each region sequentially
        for start, end, num_bins, region_name in regions:
            # Use numpy.linspace to generate boundary points for all intervals
            # N intervals require N+1 boundary points
            boundaries = np.linspace(start, end, num_bins + 1)

            # Create each sub-interval based on boundary points
            for i in range(num_bins):
                segment_start = boundaries[i]
                segment_end = boundaries[i+1]
                
                # Combine region name with interval index (starting from 1)
                indexed_region_name = f"{region_name}_{i + 1}"
                
                all_bins.append({
                    'chr': chr_name,
                    'segment_start': segment_start,
                    'segment_end': segment_end,
                    'region': indexed_region_name  # Use the new region name with index
                })

    # Convert the list to a DataFrame
    result_df = pd.DataFrame(all_bins)
    
    # Convert coordinates to integers, a standard practice in genomics
    result_df['segment_start'] = result_df['segment_start'].round().astype(int)
    result_df['segment_end'] = result_df['segment_end'].round().astype(int)

    return result_df


size_cen_df = merge_size_cen(args.size, args.cen)
segment_df = divide_chromosome_into_bins(size_cen_df)
prefix = args.cen.split('.')[0]
segment_file = f"{prefix}.segment"
segment_df.to_csv(segment_file, index=False, sep = '\t')

# Count the number of different LTR types within each interval
command_count = f"python 01.count.py --segment {segment_file} --ltr {args.ltr} --output {prefix}.count"
os.system(command_count)

# Identify the top 5 most abundant LTR types in the centromere region and aggregate others into 'others'
command_most5 = f"python 02.count_most5.py --input {prefix}.count --output {prefix}.count.most5"
os.system(command_most5)

# Merge identical intervals across different chromosomes
command_merge = f"python 03.merge.py --input {prefix}.count.most5 --output {prefix}.most5.merged"
os.system(command_merge)

# Plot line charts
command_plot = f"python 04.plot_line.py --input {prefix}.most5.merged --output {prefix}.ltr"
os.system(command_plot)