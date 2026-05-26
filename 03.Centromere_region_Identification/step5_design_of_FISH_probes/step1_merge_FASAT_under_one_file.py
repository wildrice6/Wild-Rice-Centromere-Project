import os
import argparse

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Merge multiple sequences into a single entry")
parser.add_argument('--fa', type=str, required=True, help='Path to the input FASTA file')
args = parser.parse_args()

def merge_fa_in_onefile(fa):
    """
    Merges all sequences within a single FASTA file into one continuous sequence record.
    """
    prefix = fa.rsplit('.', 1)[0]
    
    # Extract sequence data by excluding header lines
    command1 = f"grep -v '>' {fa} > temp1.fa"
    # Remove newline characters to concatenate sequences
    command2 = f"tr -d '\n' < temp1.fa > temp2.fa"
    
    os.system(command1)
    os.system(command2)
    
    # Create a new header for the merged sequence
    with open('temp3.fa', 'w') as file:
        file.write(f'>{prefix}\n')
        
    # Concatenate the header and the merged sequence data
    command3 = f"cat temp3.fa temp2.fa > {prefix}.merge.fa"
    os.system(command3)
    
    print(f'***** Successfully merged all sequences in the file into a single entry *****')
    
    # Clean up temporary files
    os.remove('temp1.fa')
    os.remove('temp2.fa')
    os.remove('temp3.fa')

if __name__ == "__main__":
    merge_fa_in_onefile(args.fa)