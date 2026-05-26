import pandas as pd
import argparse
import sys

def merge_bed_file(input_file, output_file):
    """
    Reads, sorts, and merges overlapping intervals in a BED file while considering 
    strand orientation.

    Args:
        input_file (str): Path to the input BED file.
        output_file (str): Path to the output merged BED file.
    """
    # Define column names for the BED file
    bed_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']

    try:
        print(f"Reading BED file: {input_file} ...")
        # Read the tab-separated BED file
        df = pd.read_csv(
            input_file, 
            sep='\t', 
            header=None, 
            names=bed_columns,
            # Ensure start and end coordinates are read as integers
            dtype={'chrom': str, 'start': int, 'end': int, 'name': str, 'score': float, 'strand': str}
        )
        print(f"-> Read complete; found {len(df)} records.")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
        sys.exit(1)

    # --- Step 1: Sorting ---
    # Data must be sorted by 'chrom', then 'strand', and finally 'start' for proper merging
    print("Sorting data by 'chrom', 'strand', and 'start'...")
    df_sorted = df.sort_values(by=['chrom', 'strand', 'start'])

    # --- Step 2: Merging Overlapping Intervals ---
    print("Merging overlapping intervals...")
    merged_records = []
    
    # Process the data by grouping 'chrom' and 'strand'
    for (chrom, strand), group in df_sorted.groupby(['chrom', 'strand']):
        
        # Initialize the first interval in the current group as the starting point for merging
        if group.empty:
            continue
        
        # Use .iloc[0] to retrieve the first row of data
        current_start = group.iloc[0]['start']
        current_end = group.iloc[0]['end']
        current_names = [group.iloc[0]['name']]
        current_scores = [group.iloc[0]['score']]

        # Iterate through intervals starting from the second record
        for i in range(1, len(group)):
            next_row = group.iloc[i]
            next_start = next_row['start']
            
            # Check for overlap (start of the next interval occurs before the end of the current merged interval)
            if next_start < current_end:
                # Overlap detected; extend the current merged interval endpoint if necessary
                current_end = max(current_end, next_row['end'])
                current_names.append(next_row['name'])
                current_scores.append(next_row['score'])
            else:
                # No overlap; finalize the previous merged interval and append it to the results
                merged_records.append({
                    'chrom': chrom,
                    'start': current_start,
                    'end': current_end,
                    'name': ','.join(map(str, sorted(list(set(current_names))))), # Deduplicate and join names
                    'score': max(current_scores), # Retrieve the maximum bitscore
                    'strand': strand
                })
                
                # Initialize a new merged interval
                current_start = next_row['start']
                current_end = next_row['end']
                current_names = [next_row['name']]
                current_scores = [next_row['score']]

        # Append the final merged interval of the group to the results
        merged_records.append({
            'chrom': chrom,
            'start': current_start,
            'end': current_end,
            'name': ','.join(map(str, sorted(list(set(current_names))))),
            'score': max(current_scores),
            'strand': strand
        })

    # --- Step 3: Saving Results ---
    if not merged_records:
        print("No records available to merge or output.")
        # Create an empty file
        open(output_file, 'w').close()
        return
        
    print(f"-> Merging complete; generated {len(merged_records)} new records.")
    # Convert the results list to a DataFrame
    df_merged = pd.DataFrame(merged_records)
    
    # Re-sort the final results (optional, but recommended for consistency)
    df_merged_sorted = df_merged.sort_values(by=['chrom', 'start'])
    
    print(f"Saving results to: {output_file} ...")
    df_merged_sorted.to_csv(output_file, sep='\t', header=False, index=False)
    
    print("\nProcessing successful!")

def main():
    """
    Main function to parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Sorts and merges overlapping intervals in a BED file while accounting for strand orientation.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', "--input_file", 
        help="Input BED filename.\nExample: input.bed"
    )
    parser.add_argument(
        '-o', "--output_file", 
        help="Output merged BED filename.\nExample: merged.bed"
    )
    
    args = parser.parse_args()
    merge_bed_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()