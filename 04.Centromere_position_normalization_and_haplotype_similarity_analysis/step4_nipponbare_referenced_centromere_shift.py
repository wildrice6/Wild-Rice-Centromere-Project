import pandas as pd
import argparse
import sys

def main():
    # 1. Argument parsing setup
    parser = argparse.ArgumentParser(description="Calculate the positional differences between each assembly and the reference AA_Osat_hap2")
    parser.add_argument("--input", required=True, help="Path to the input TSV file")
    parser.add_argument("--output", required=True, help="Path to the output TSV result file")
    
    args = parser.parse_args()

    try:
        # 2. Load data
        df = pd.read_csv(args.input, sep='\t')
        
        # Check for required columns
        required_cols = ['ChrID', 'Normalized_Start', 'Normalized_End', 'assembly']
        for col in required_cols:
            if col not in df.columns:
                print(f"Error: Missing required column '{col}' in the input file")
                sys.exit(1)

        # 3. Calculate the midpoint for each entry
        df['Midpoint'] = (df['Normalized_Start'] + df['Normalized_End']) / 2

        # 4. Extract reference data (AA_Osat_hap2)
        ref_assembly = "AA_Osat_hap2"
        ref_df = df[df['assembly'] == ref_assembly].copy()

        if ref_df.empty:
            print(f"Error: Reference assembly '{ref_assembly}' not found in the input file")
            sys.exit(1)

        # Create a reference dictionary for efficient lookup: {ChrID: Midpoint}
        # If a chromosome has multiple entries, take the first one (genome data typically has one interval per chromosome)
        ref_map = ref_df.groupby('ChrID')['Midpoint'].first().to_dict()

        # 5. Calculate differences
        results = []
        
        for index, row in df.iterrows():
            chr_id = row['ChrID']
            current_mid = row['Midpoint']
            vs_assembly = row['assembly']
            
            if chr_id in ref_map:
                diff = current_mid - ref_map[chr_id]
                results.append({
                    'ChrID': chr_id,
                    'difference': diff,
                    'VS_assembly': vs_assembly
                })
            else:
                # Skip or report error if the chromosome is missing in the reference assembly
                print(f"Warning: Chromosome {chr_id} not found in reference {ref_assembly}, skipping.")

        # 6. Convert to DataFrame and save
        output_df = pd.DataFrame(results)
        
        # Sort by VS_assembly and ChrID for better readability
        output_df = output_df.sort_values(by=['VS_assembly', 'ChrID'])

        # Save as TSV
        output_df.to_csv(args.output, sep='\t', index=False)
        print(f"Calculation complete! Results saved to: {args.output}")

    except Exception as e:
        print(f"An error occurred during processing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()