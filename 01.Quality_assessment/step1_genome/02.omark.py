import os
import argparse
import pandas as pd
parser = argparse.ArgumentParser(description="omark评估")
parser.add_argument('--faa', type=str, required=True, help='蛋白质序列文件')
args = parser.parse_args()

def omark(faa):
    prefix = faa.rsplit('.', 1)[0]
    command1 = f"omamer search --db ~/mysoftware/LUCA.h5 --query {faa} --out {prefix}.omamer"
    os.system(command1)
    command2 = f"omark -f {prefix}.omamer -d ~/mysoftware/LUCA.h5 -o {prefix}_output"
    os.system(command2)
    print(f'已完成{faa}的omark评估')

omark(args.faa)