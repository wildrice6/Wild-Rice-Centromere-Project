import os
import argparse
import pandas as pd
parser = argparse.ArgumentParser(description="Map sequencing reads to the reference genome")
parser.add_argument('-d', '--data', type=str, required=True, help='Sequencing data type: ONT or HiFi')
parser.add_argument('-s', '--seq', type=str, required=True, help='Sequencing data type: ONT or HiFi')
parser.add_argument('-g', '--genome', type=str, required=True, help='Genomic data')
args = parser.parse_args()

def hifi_mapping_coverage(hifi, genome):
    prefix = hifi.split('.')[0]
    command1 = f"bedtools bamtofastq -i {hifi} -fq {prefix}.hifi.fastq"
    os.system(command1)
    command2 = f"minimap2 -ax map-pb {genome} {prefix}.hifi.fastq -t 64 > {prefix}.hifi.sam"
    os.system(command2)
    command3 = f"samtools view -@ 64 -Sb {prefix}.hifi.sam > {prefix}.hifi.bam"
    os.system(command3)
    command4 = f"samtools sort -@ 64 -o {prefix}.hifi.sorted.bam {prefix}.hifi.bam"
    os.system(command4)

def ont_mapping_coverage(ont, genome):
    prefix = ont.split('.')[0]
    command1 = f"minimap2 --secondary=no -ax map-ont {genome} {ont} -t 64 > {prefix}.ont.sam"
    command2 = f"samtools view -@ 64 -Sb {prefix}.ont.sam > {prefix}.ont.bam"
    os.system(command2)
    command3 = f"samtools sort -@ 64 -o {prefix}.ont.sorted.bam {prefix}.ont.bam"
    os.system(command3)


if args.data == 'hifi':
    hifi_mapping_coverage(args.seq, args.genome)
elif args.data == 'ont':
    ont_mapping_coverage(args.seq, args.genome)