import os
import numpy as np
import pandas as pd
import sys
sys.path.append('/share/org/YZWL/yzwl_shahd/mysoftware/frompip')
import miniseq
import argparse

# Parse command-line arguments
argparser = argparse.ArgumentParser(description='Obtain the consensus sequence from input data')
argparser.add_argument('-i', '--input', type=str, required=True)
args = argparser.parse_args()
inputfile = args.input
fileprefix = inputfile.rsplit('.', 1)[0]

satellites = miniseq.FASTA(filename=inputfile)

# Load the input file
seqdata = []
for satellite in satellites:
    seqdata.append(satellite.get_sequence())

seq_list = [list(string) for string in seqdata]
seq_df = pd.DataFrame(seq_list)

# Generate consensus sequence by identifying the most frequent base at each position
consensus_seq = seq_df.apply(lambda x: x.value_counts().idxmax())
consensus_str = consensus_seq.astype(str)
conbined_str = consensus_seq.str.cat(sep='')

outputfile = fileprefix + '.consensus'
seq_identifier = '>' + fileprefix + '\n'

conbined_str = seq_identifier + conbined_str
with open(outputfile, 'w') as file:
    file.write(f"{conbined_str}\n")

print(f'Consensus sequence for {fileprefix} has been successfully extracted.')
print('*'*50)