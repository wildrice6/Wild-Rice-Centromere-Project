import os
import pandas as pd
import argparse

# Note: Prior to execution, the gff2bed utility from the BEDOPS software suite 
# must be used via the command line to convert the original GFF output into a BED file.

parser = argparse.ArgumentParser(description="Standardize a BED file derived from GFF conversion")
parser.add_argument('--bed', type=str, required=True, help='Path to the input BED file')
args = parser.parse_args()
bed = args.bed

def get_standard_bed(bed):
    """
    Input a BED file converted from GFF format and generate a standardized 6-column BED file.
    """
    material = bed.split('.')[0]
    
    # Read the file using both tab and semicolon as delimiters to handle GFF-derived attributes
    temp = pd.read_csv(bed, header=None, sep='\t|;', engine='python')
    
    # Select specific columns to form a standard BED structure: 
    # Chromosome (0), Start (1), End (2), Name (9), Score (4), and Strand (5)
    output = temp.loc[:, [0, 1, 2, 9, 4, 5]]
    
    # Clean the 'Name' column by removing the 'Name=' prefix
    output[9] = output[9].replace(r'Name=', '', regex=True)
    
    # Optional formatting (currently commented out in original implementation):
    # output[9] = output[9].replace(r'rRNA', 'rDNA', regex=True)
    # output[9] = output[9] + '_' + (output[9].index + 1).astype(str)
    # output[9] = material + '-' + output[9]
    
    output_name = f"{bed}.txt"
    output.to_csv(output_name, header=False, sep='\t', index=False)
    
    print(f"***** File {bed} has been successfully converted to a standardized BED format *****")

if __name__ == "__main__":
    get_standard_bed(bed)