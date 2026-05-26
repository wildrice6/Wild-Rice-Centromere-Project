import os
import argparse
import pandas as pd

# Initialize argument parser for rDNA annotation
parser = argparse.ArgumentParser(description="Perform rDNA annotation on genomic sequences")
parser.add_argument('-fa', '--fasta', type=str, required=True, help='Path to the genome FASTA file')
args = parser.parse_args()

def barrnap(fa):
    """
    Input a FASTA sequence file and output the rDNA annotation in GFF format.
    """
    # Extract the file prefix by removing the extension
    prefix = fa.rsplit('.', 1)[0]
    
    # Construct the command for pybarrnap using the eukaryotic kingdom and accurate mode
    command1 = f"pybarrnap -k euk {fa} --accurate > {prefix}.rDNA.gff"
    
    # Execute the command
    os.system(command1)
    
    # Output completion status
    print(f"***** rDNA annotation for {fa} has been completed *****")

if __name__ == "__main__":
    barrnap(args.fasta)