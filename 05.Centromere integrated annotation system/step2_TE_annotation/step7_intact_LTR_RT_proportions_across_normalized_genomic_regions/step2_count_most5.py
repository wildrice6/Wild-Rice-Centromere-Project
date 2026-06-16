#!/usr/bin/env python3
"""
TSV Data Processing Script: Handling LTR Transposon Data.
Retains the top 5 most abundant non-unknown LTR types within the centromere region, 
aggregating all other types into an "others" category.
"""

import pandas as pd
import argparse
import sys
from collections import Counter


def identify_top_ltr_types(df, top_n=5):
    """
    Identifies the top N LTR types by abundance within the centromere region (excluding "unknown" types).
    
    Args:
        df: DataFrame
        top_n: Number of LTR types to retain
    
    Returns:
        list: Column names of the top N LTR types
    """
    # Filter for centromere region
    centromere_mask = df['region'].str.contains('centromere_', na=False)
    centromere_df = df[centromere_mask]
    
    if centromere_df.empty:
        print("Warning: No data found for the centromere region")
        return []
    
    # Get all LTR column names (excluding types ending with "_unknown")
    ltr_columns = [col for col in df.columns 
                   if col.startswith('LTR_') and not col.endswith('_unknown')]
    
    # Calculate total abundance of each LTR type within the centromere region
    ltr_counts = {}
    for col in ltr_columns:
        total_count = centromere_df[col].sum()
        if total_count > 0:  # Only consider LTR types with data
            ltr_counts[col] = total_count
    
    # Sort by abundance and retrieve the top N types
    top_ltr_types = sorted(ltr_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
    top_ltr_columns = [item[0] for item in top_ltr_types]
    
    print(f"Top {top_n} LTR types in the centromere region:")
    for col, count in top_ltr_types:
        print(f"  {col}: {count}")
    
    return top_ltr_columns


def process_tsv(input_file, output_file, top_n=5):
    """
    Processes the TSV file.
    
    Args:
        input_file: Path to the input file
        output_file: Path to the output file
        top_n: Number of LTR types to retain
    """
    try:
        # Read the TSV file
        df = pd.read_csv(input_file, sep='\t')
        print(f"Successfully read file: {input_file}")
        print(f"Data dimensions: {df.shape}")
        
        # Verify required columns
        required_columns = ['chr', 'segment_start', 'segment_end', 'region']
        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f"Missing required column: {col}")
        
        # Get all LTR column names
        all_ltr_columns = [col for col in df.columns if col.startswith('LTR_')]
        print(f"Found {len(all_ltr_columns)} LTR columns")
        
        # Identify the top N LTR types within the centromere region
        top_ltr_columns = identify_top_ltr_types(df, top_n)
        
        if not top_ltr_columns:
            print("Warning: No valid LTR types found; aggregating all LTRs into 'others'")
        
        # Create a new DataFrame
        result_df = df[['chr', 'segment_start', 'segment_end', 'region']].copy()
        
        # Add columns for the top N LTR types
        for col in top_ltr_columns:
            result_df[col] = df[col]
        
        # Calculate 'others' column: sum of all other LTR types
        other_ltr_columns = [col for col in all_ltr_columns if col not in top_ltr_columns]
        
        if other_ltr_columns:
            result_df['others'] = df[other_ltr_columns].sum(axis=1)
            print(f"Aggregating {len(other_ltr_columns)} LTR types into 'others':")
            for col in other_ltr_columns:
                print(f"  {col}")
        else:
            result_df['others'] = 0
            print("No other LTR types to aggregate into 'others'")
        
        # Save results
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"Processing complete; results saved to: {output_file}")
        print(f"Output file dimensions: {result_df.shape}")
        
        # Display output file column names
        print("Columns in the output file:")
        for col in result_df.columns:
            print(f"  {col}")
            
    except FileNotFoundError:
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during processing: {e}")
        sys.exit(1)


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Process LTR transposon TSV data, retaining the top 5 most abundant LTR types in the centromere region'
    )
    
    parser.add_argument(
        '--input', 
        required=True, 
        help='Path to the input TSV file'
    )
    
    parser.add_argument(
        '--output', 
        required=True, 
        help='Path to the output TSV file'
    )
    
    parser.add_argument(
        '--top_n', 
        type=int, 
        default=5, 
        help='Number of LTR types to retain (default: 5)'
    )
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("LTR Transposon Data Processing Script")
    print("=" * 60)
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Retained LTR types: {args.top_n}")
    print("-" * 60)
    
    process_tsv(args.input, args.output, args.top_n)


if __name__ == "__main__":
    main()