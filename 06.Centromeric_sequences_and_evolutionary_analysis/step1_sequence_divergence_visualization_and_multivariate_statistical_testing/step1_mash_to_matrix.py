import pandas as pd
import argparse
import os
import sys

def get_clean_name(path):
    """
    Extract the core identifier (e.g., Species_Chromosome) from the file path.
    Example: 01.mash/AA_Ogla_hap2_Chr03.cenRegion.fasta -> AA_Ogla_hap2_Chr03
    """
    # 1. Extract the filename (remove directory path)
    filename = os.path.basename(path)
    
    # 2. Strip extensions
    # Specifically target the '.cenRegion.fasta' suffix or standard genomic extensions
    if ".cenRegion.fasta" in filename:
        clean_name = filename.replace(".cenRegion.fasta", "")
    elif ".fasta" in filename:
        clean_name = filename.replace(".fasta", "")
    elif ".fa" in filename:
        clean_name = filename.replace(".fa", "")
    else:
        # Fallback: take the substring before the first dot if extensions do not match
        clean_name = filename.split('.')[0]
        
    return clean_name

def main():
    parser = argparse.ArgumentParser(description="Convert 'mash dist' output into a square distance matrix.")
    parser.add_argument("input_file", help="Path to the mash dist output file (TSV format).")
    parser.add_argument("output_file", help="Path to save the generated matrix (TSV format).")
    
    args = parser.parse_args()
    
    input_path = args.input_file
    output_path = args.output_file

    print(f"Reading input file: {input_path} ...")
    
    try:
        # 1. Read Mash dist output
        # Mash output typically lacks headers; column order is defined manually
        # Columns: Reference, Query, Distance, P-value, Matching-hashes
        df = pd.read_csv(input_path, sep='\t', header=None, names=['Ref', 'Query', 'Dist', 'Pval', 'Hashes'])
        
        # 2. Sanitize names
        # Apply the cleaning function to both 'Ref' and 'Query' columns
        df['Ref_Clean'] = df['Ref'].apply(get_clean_name)
        df['Query_Clean'] = df['Query'].apply(get_clean_name)
        
        print("Transforming to matrix format...")
        
        # 3. Pivot the table into a square matrix
        # Index = Row names, Columns = Column names, Values = Distance scores
        matrix = df.pivot(index='Ref_Clean', columns='Query_Clean', values='Dist')
        
        # 4. Handle missing values
        # Mash usually provides full pairwise comparisons; missing values (e.g., diagonals) are filled with 0
        matrix.fillna(0, inplace=True)
        
        # 5. Save results
        # Export as a tab-separated file including row and column headers
        matrix.to_csv(output_path, sep='\t')
        
        print(f"Success! Matrix saved to: {output_path}")
        print(f"Matrix shape: {matrix.shape[0]} rows x {matrix.shape[1]} columns")
        
    except Exception as e:
        print(f"\n[Error] An error occurred: {e}")
        print("Ensure the input is a valid 'mash dist' output and that the pandas library is installed.")
        sys.exit(1)

if __name__ == "__main__":
    main()