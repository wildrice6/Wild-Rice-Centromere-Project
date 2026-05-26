#!/usr/bin/env python3
"""
Batch Merging of Segmentation Results
=====================================
Merge TSV files from different chromosomes belonging to the same hap & rep into a single file.
Input file naming convention:
    FF_Obra_Chr01.hap1.rep2.segmentation.tsv
Output file example:
    FF_Obra_hap1.rep2.segmentation.merged.tsv

Usage:
    python merge_segmentation.py        # Execute directly within the directory containing the TSV files.
"""

import glob, os, re, collections
import pandas as pd

# -------- 1. Collect all segmentation files -----------------------------
SEG_FILES = glob.glob("FF_Obra_Chr*.hap*.rep*.segmentation.tsv")
if not SEG_FILES:
    raise SystemExit("❌ No segmentation TSV files found in the directory. Please check the path or naming format.")

# Use (hap, rep) as keys to collect the corresponding file lists into a dictionary
pat = re.compile(r"FF_Obra_Chr(\d+)\.hap(\d+)\.rep(\d+)\.segmentation\.tsv$")
grp = collections.defaultdict(list)

for f in SEG_FILES:
    m = pat.search(f)
    if not m:          # Skip files that do not match the naming convention
        continue
    chr_id, hap, rep = m.groups()
    grp[(hap, rep)].append((int(chr_id), f))   # Record chromosome index for sorting

# -------- 2. Merge each (hap, rep) combination -------------------------------
for (hap, rep), lst in grp.items():
    # Sort by chromosome number in ascending order
    lst.sort()        # Elements in lst are (chr_num, file)
    dfs = []
    for chr_num, f in lst:
        df = pd.read_csv(f, sep="\t")
        if "chrom" not in df.columns:          # Add 'chrom' column if it is missing
            df.insert(0, "chrom", f"Chr{chr_num:02d}")
        dfs.append(df)

    merged = pd.concat(dfs, ignore_index=True)
    out_name = f"FF_Obra_hap{hap}.rep{rep}.segmentation.merged.tsv"
    merged.to_csv(out_name, sep="\t", index=False)

    print(f"✅ {hap=}, {rep=} → Merged {len(lst)} files, wrote {out_name} "
          f"({len(merged)} lines)")

print(f"🎉 All tasks completed, generated {len(grp)} merged files in total.")