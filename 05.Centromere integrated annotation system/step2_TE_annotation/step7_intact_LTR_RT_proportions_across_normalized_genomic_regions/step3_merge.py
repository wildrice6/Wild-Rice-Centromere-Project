#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse

def main():
    """
    Main function: Read data, aggregate LTR counts by 'region', 
    calculate the proportion of each LTR type within its respective region, 
    and save the results.
    """
    # 1. Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description="Summarize LTR counts by region and convert them to proportions.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to the input TSV file with LTR counts per segment.'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path for the output TSV file with LTR proportions.'
    )

    args = parser.parse_args()

    # 2. Read the input file
    try:
        print(f"--- Reading data from: {args.input} ---")
        df = pd.read_csv(args.input, sep='\t')
        
        # Correct potential column name spelling errors (e.g., 'hr' -> 'chr')
        if df.columns[0].lower() == 'hr':
            df.rename(columns={df.columns[0]: 'chr'}, inplace=True)

        print("Successfully loaded the data. Preview:")
        print(df.head())

    except FileNotFoundError:
        print(f"Error: Input file not found at '{args.input}'")
        return
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return

    # 3. Step 1: Aggregate LTR counts by 'region'
    print("\n--- Step 1: Aggregating LTR counts by 'region' ---")

    if 'region' not in df.columns:
        print("Error: The input file must contain a 'region' column.")
        return

    # Identify all LTR count columns (starting from the 5th column)
    count_columns = df.columns[4:]

    # Group by 'region' and sum all count columns
    aggregated_counts_df = df.groupby('region')[list(count_columns)].sum()

    print("Aggregation complete. Preview of total counts:")
    print(aggregated_counts_df.head())

    # 4. Step 2: Convert counts to proportions
    print("\n--- Step 2: Converting counts to proportions within each region ---")
    
    # Calculate the total number of LTRs for each row (each region)
    # axis=1 indicates summation across rows
    region_totals = aggregated_counts_df.sum(axis=1)

    # Use the .div() method to divide each row value by the row total
    # axis=0 ensures that the index of region_totals (a Series) aligns with the row index of aggregated_counts_df
    proportions_df = aggregated_counts_df.div(region_totals, axis=0)
    
    # Handle division by zero: If a region's total LTR count is 0, the division results in NaN (Not a Number).
    # Fill these NaN values with 0 since the actual proportion is 0.
    proportions_df.fillna(0, inplace=True)
    
    # Reset 'region' from index to a standard column
    proportions_df.reset_index(inplace=True)
    
    print("Proportion calculation complete. Preview of the final data:")
    print(proportions_df.head())

    # 5. Save results to the output file
    try:
        print(f"\n--- Saving proportional data to: {args.output} ---")
        # Specify floating-point format for better readability
        proportions_df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
        print("Process completed successfully!")

    except Exception as e:
        print(f"An error occurred while saving the file: {e}")

if __name__ == "__main__":
    main()