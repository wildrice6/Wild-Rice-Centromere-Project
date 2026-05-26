import os
import argparse

argparser = argparse.ArgumentParser(description='Convert tabular data into sequence files')
argparser.add_argument('-i', '--input', type=str, required=True)
argparser.add_argument('-m', '--model', type=str, default='nucl')
args = argparser.parse_args()
fasta = args.input
model = args.model

pwd = os.getcwd()
name = pwd.split('/')[-1]
prefix = fasta.rsplit('.', 1)[0]
shname = f"{prefix}.sh"

command1 = f'#!/bin/bash\n#CSUB -J {name}\n#CSUB -q c01\n#CSUB -o {prefix}.out\n#CSUB -e test.out\n#CSUB -n 4\n#CSUB -cwd {pwd}\n\nsource /share/apps/anaconda3/bin/activate\nconda activate ga'
command2 = f"\nmafft --auto --thread 4 --anysymbol --maxiterate 1000 {fasta} > {prefix}.mafft\n"
command3 = f"trimal -in {prefix}.mafft -out {prefix}.trimal -automated1\n"
command4_nucl = f"fasttree -nt -gtr {prefix}.trimal > {prefix}.tree\n"
command4_prot = f"fasttree {prefix}.trimal > {prefix}.tree\n"
command5 = f"\nmkdir {prefix}\n"
command6 = f"\nmv {prefix}.mafft {prefix}.trimal {prefix}.tree {prefix}\n"

if model == 'nucl':
        with open(shname, 'w') as file:
                file.write(command1)
                file.write(command2)
                file.write(command3)
                file.write(command4_nucl)
                file.write(command5)
                file.write(command6)
elif model == 'prot':
        with open(shname, 'w') as file:
                file.write(command1)
                file.write(command2)
                file.write(command3)
                file.write(command4_prot)
                file.write(command5)
                file.write(command6)