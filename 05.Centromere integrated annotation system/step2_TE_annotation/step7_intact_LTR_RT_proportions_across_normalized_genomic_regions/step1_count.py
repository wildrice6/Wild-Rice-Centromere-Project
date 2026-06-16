#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
from tqdm import tqdm

def main():
    """
    Main function: Parse arguments, perform LTR counting, and save the results.
    """
    # 1. Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description="Count different types of LTRs within genomic segments.",
        formatter_class=argparse.RawTextHelpFormatter # Maintain help message formatting
    )
    parser.add_argument('--segment', required=True, help='Path to the segment definition file (TSV format).')
    parser.add_argument('--ltr', required=True, help='Path to the LTR annotation file (TSV/space-separated format).')
    parser.add_argument('--output', required=True, help='Path for the output TSV file.')
    
    args = parser.parse_args()

    # 2. Read input files
    print("--- Reading input files ---")
    try:
        # Read the segment file, assuming columns are 'chr', 'segment_start', 'segment_end', 'region'
        df_segment = pd.read_csv(args.segment, sep='\t')
        print(f"Successfully loaded {len(df_segment)} segments from {args.segment}")

        # Read the LTR file; based on typical annotation formats, it lacks a header and may be space-delimited
        # Focus is placed on the first four columns: chromosome, start, end, and type
        ltr_cols = ['chr', 'start', 'end', 'type', 'id', 'strand', 'family']
        df_ltr = pd.read_csv(args.ltr, sep='\s+', header=None, names=ltr_cols, usecols=[0, 1, 2, 3])
        print(f"Successfully loaded {len(df_ltr)} LTRs from {args.ltr}")

    except FileNotFoundError as e:
        print(f"Error: Input file not found - {e}")
        return
    except Exception as e:
        print(f"An error occurred while reading files: {e}")
        return

    # 3. Core Logic: Count LTRs within each segment
    print("\n--- Counting LTRs in each segment ---")
    
    # To optimize calculation, pre-calculate the midpoint of each LTR.
    # Assigning an LTR to a segment based on its midpoint is a robust method
    # that prevents a single LTR from being redundantly counted in adjacent segments.
    df_ltr['midpoint'] = (df_ltr['start'] + df_ltr['end']) // 2

    # Identify all unique LTR types to be used as column names in the output
    all_ltr_types = sorted(df_ltr['type'].unique())
    print(f"Found {len(all_ltr_types)} unique LTR types.")

    # Initialize an empty DataFrame to store counts, with LTR types as columns
    count_matrix = pd.DataFrame(0, index=df_segment.index, columns=all_ltr_types)

    # Iterate through each segment to perform statistics
    # tqdm is utilized as a progress bar for monitoring large-scale data processing
    # If not installed, run: pip install tqdm
    for index, row in tqdm(df_segment.iterrows(), total=len(df_segment), desc="Processing Segments"):
        chr_name = row['chr']
        seg_start = row['segment_start']
        seg_end = row['segment_end']

        # Filter LTRs located on the same chromosome whose midpoints fall within the segment interval
        # Use a standard half-open interval [seg_start, seg_end)
        mask = (df_ltr['chr'] == chr_name) & \
               (df_ltr['midpoint'] >= seg_start) & \
               (df_ltr['midpoint'] < seg_end)
        
        ltrs_in_segment = df_ltr.loc[mask]

        # If LTRs are present in the segment
        if not ltrs_in_segment.empty:
            # Count the occurrences of each unique LTR type
            counts = ltrs_in_segment['type'].value_counts()
            
            # Update the counting results in the matrix
            count_matrix.loc[index, counts.index] = counts.values

    # 4. Merge original segment data with the calculated counting results
    print("\n--- Merging results ---")
    df_final = pd.concat([df_segment, count_matrix], axis=1)

    # 5. Save results to the output file
    try:
        df_final.to_csv(args.output, sep='\t', index=False)
        print(f"\nSuccessfully saved the result to {args.output}")
        print("--- Final DataFrame Head ---")
        print(df_final.head())
    except Exception as e:
        print(f"An error occurred while saving the file: {e}")


if __name__ == "__main__":
    main()