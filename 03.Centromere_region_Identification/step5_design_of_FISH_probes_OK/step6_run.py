import os
import argparse
import pandas as pd

# Define command-line arguments for the probe design pipeline
parser = argparse.ArgumentParser(description="Probe Design Pipeline")
parser.add_argument('--target_fa', type=str, required=True, help='Target FASTA sequence file')
parser.add_argument('--genome', type=str, required=True, help='Reference genome sequence file')
parser.add_argument('--cen', type=str, required=True, help='Centromere regions coordinate file')
args = parser.parse_args()

def main(target_fa, genome, cen):
    # Merge all sequences in the FASTA file into a single entry to facilitate k-mer partitioning
    command_merge = f"python ~/quickstart/merge_FASAT_under_one_file.py --fa {target_fa}"
    os.system(command_merge)

    # Partition the target FASTA sequence into k-mers
    target_fa_prefix = target_fa.rsplit('.', 1)[0]
    command_kmer = f"python get_kmer.py --input {target_fa_prefix}.merge.fa --output {target_fa_prefix}.merge.kmer"
    os.system(command_kmer)

    # Split the multi-sequence k-mer file into individual FASTA files (one sequence per file)
    command_split_kmer = f"python split_kmer.py --input {target_fa_prefix}.merge.kmer --output {target_fa_prefix}.kmer"
    os.system(command_split_kmer)

    # Align the partitioned k-mers back to the reference genome using BLASTN
    # Output files will be located in the directory: {genome_prefix}.{target_fa_prefix}.kmer
    genome_prefix = genome.rsplit('.')[0]
    command_blastn = f"for file in {target_fa_prefix}.kmer/*.fasta; do python ~/quickstart/blast.py -m short_blastn -q ${{file}} -t {genome} ; done"
    os.system(command_blastn)

    # Visualize the distribution of probes across the genome
    # First, index the genome using samtools
    command_fai = f"samtools faidx {genome}"
    os.system(command_fai)
    
    # Generate distribution plots for each alignment result
    command_plot = f"for file in {genome_prefix}.{target_fa_prefix}.kmer/*.txt; do python probe_plot.py --blast_output ${{file}} --genome_size {genome}.fai --cen {cen} --output ${{file}}.png; done"
    os.system(command_plot)

if __name__ == "__main__":
    main(args.target_fa, args.genome, args.cen)