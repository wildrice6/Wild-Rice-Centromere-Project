import pandas as pd
import argparse
import os

def get_material_name(filepath):
    """
    Extracts the material name from the file path.
    New rule: The first substring separated by a period.
    Example: Extracts 'SampleA' from 'SampleA.some.data.stat'.
    """
    # 1. Retrieve the filename excluding the path
    filename = os.path.basename(filepath)
    # 2. Split by the period '.' and retrieve the first component
    material_name = filename.split('.')[0]
    return material_name

def main():
    """
    Main function for parameter parsing and file processing.
    """
    # --- 1. Set up command-line argument parsing ---
    parser = argparse.ArgumentParser(description="Aggregates block_size for headerless TSV files, calculates percentages, and appends a material name column.")
    parser.add_argument("--input", required=True, help="Path to the input TSV file (no header).")
    parser.add_argument("--output", required=True, help="Path to the output TSV file.")
    
    args = parser.parse_args()

    # --- 2. Check if the input file exists ---
    if not os.path.exists(args.input):
        print(f"Error: Input file does not exist -> {args.input}")
        return

    try:
        # --- 3. Extract material name (using updated logic) ---
        material_name = get_material_name(args.input)
        print(f"Material name extracted from filename '{os.path.basename(args.input)}': {material_name}")

        # --- 4. Read the TSV file using pandas (no header) ---
        df = pd.read_csv(
            args.input, 
            sep='\s+', 
            header=None, 
            names=['block_size', 'count']
        )

        # --- 5. Group by 'block_size' and aggregate 'count' via summation ---
        print("Potential duplicate block_size entries detected; performing aggregation...")
        aggregated_df = df.groupby('block_size')['count'].sum().reset_index()
        
        # --- 6. Calculate the total sum of the aggregated 'count' column ---
        total_count = aggregated_df['count'].sum()
        print(f"Total aggregated 'count' value: {total_count}")

        # --- 7. Calculate percentage and create a new column ---
        if total_count > 0:
            aggregated_df['percentage'] = (aggregated_df['count'] / total_count) * 100
        else:
            aggregated_df['percentage'] = 0.0

        # --- 8. Add material name column ---
        aggregated_df['material_name'] = material_name

        # --- 9. Adjust column order ---
        final_df = aggregated_df[['material_name', 'block_size', 'count', 'percentage']]

        # --- 10. Save the results to the output file ---
        final_df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
        
        print(f"Processing complete. Results saved to: {args.output}")

    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == "__main__":
    main()