import os
import argparse
import pandas as pd

# Initialize argument parser for rDNA quantification
parser = argparse.ArgumentParser(description="Calculate the counts of various rDNA types")
parser.add_argument('--bed', type=str, required=True, help='Path to the input BED file')
args = parser.parse_args()
bed = args.bed

def get_count(bed):
    """
    Quantify the occurrences of different rDNA features from the input BED file.
    """
    prefix = bed.rsplit('.', 1)[0]
    species = bed.split('.')[0]
    
    # Read the standardized BED file
    temp = pd.read_csv(bed, header=None, sep='\t')
    
    # Ensure the 'Name' prefix is removed if present (safety check)
    temp = temp.replace(r'Name=', '', regex=True)
    
    # Group by chromosome/scaffold (column 0) and count occurrences of each rDNA type (column 3)
    count = temp.groupby(0)[3].value_counts().unstack(fill_value=0)
    
    # Add species metadata
    count['species'] = species
    
    # Define output columns and export results
    output_name = f"{prefix}.num"
    # Note: Column labels correspond to standard barrnap outputs
    output = count.loc[:, ['18S_rRNA', '5_8S_rRNA', '28S_rRNA', '5S_rRNA', 'species']]
    output.to_csv(output_name, sep='\t')
    
    print(f'***** rDNA counts for {bed} have been successfully calculated *****')

if __name__ == "__main__":
    get_count(bed)