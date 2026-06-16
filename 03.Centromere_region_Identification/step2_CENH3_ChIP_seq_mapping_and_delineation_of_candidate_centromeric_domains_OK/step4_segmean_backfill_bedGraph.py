#!/usr/bin/env python3
"""
Batch Mapping of Segmentation Mean Values to Multiple bedGraph Files
--------------------------------------------------------------------
Input:
- A directory containing multiple bedGraph files (e.g., AA_Olon_hap1.sample1.CENH3.bdg)
- A directory containing multiple segmentation files (e.g., AA_Olon_hap1.rep1.segmentation.merged.tsv)

Functionality:
- Automatically identify the matching relationship between hapX.sampleY and hapX.repY (Y must correspond exactly).
- Map the 'seg.mean' from segmentation intervals to the bedGraph files.
- Generate new *.segmean.bdg files.

Usage:
python batch_map_segmean.py \
    --bdg_dir AA_Olon_bdg_files \
    --seg_dir AA_Olon_segmentation_merge \
    --out_dir AA_Olon_segmentation_bdg_files \
    --prefix AA_Olon \
    --threads 4
"""

import argparse, re, sys
from pathlib import Path
import importlib.util

# ---------------- 1. Argument Parsing ----------------
cli = argparse.ArgumentParser()
cli.add_argument("--bdg_dir", required=True, type=Path, help="Directory containing *.bdg files")
cli.add_argument("--seg_dir", required=True, type=Path, help="Directory containing segmentation files")
cli.add_argument("--out_dir", required=True, type=Path, help="Output directory (created automatically)")
cli.add_argument("--prefix", default="", help="Filename prefix (e.g., AA_Olon); leave empty for wildcard matching")
cli.add_argument("--threads", type=int, default=1, help="Number of parallel threads")
args = cli.parse_args()

args.out_dir.mkdir(parents=True, exist_ok=True)

# ---------------- 2. Dynamic Loading of Underlying Modules ----------------
seg_map_file = Path(__file__).with_name("map_segmean_to_bdg.py")
if not seg_map_file.exists():
    sys.exit(f"[ERROR] Cannot find map_segmean_to_bdg.py: {seg_map_file}")
spec = importlib.util.spec_from_file_location("map_seg", seg_map_file)
map_seg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(map_seg)

# ---------------- 3. Constructing Regular Expressions and Indices ----------------
pre = re.escape(args.prefix) + "_" if args.prefix else r".*?_"
pat_bdg = re.compile(fr"{pre}hap(\d+)\.sample(\d+)\.CENH3\.bdg$", re.I)
pat_seg = re.compile(fr"{pre}hap(\d+)\.rep(\d+)\.segmentation\.merged\.tsv$", re.I)

seg_index = {}
for seg_path in args.seg_dir.glob("*.tsv"):
    m = pat_seg.match(seg_path.name)
    if m:
        hap, rep = m.groups()
        seg_index[(hap, rep)] = seg_path

# ---------------- 4. Matching bedGraph and Segmentation Files ----------------
tasks = []
for bdg_path in args.bdg_dir.glob("*.bdg"):
    m = pat_bdg.match(bdg_path.name)
    if not m:
        print(f"[WARN] Filename does not match the naming convention, skipping: {bdg_path.name}")
        continue
    hap, sample = m.groups()
    seg_path = seg_index.get((hap, sample))
    if seg_path is None:
        print(f"[WARN] No matching segmentation file found for: hap{hap}, rep{sample}. Skipping {bdg_path.name}")
        continue
    out_name = bdg_path.with_suffix("").name + ".segmean.bdg"
    out_path = args.out_dir / out_name
    tasks.append((seg_path, bdg_path, out_path))

if not tasks:
    sys.exit("[ERROR] ❌ No matching (bdg, seg) pairs found. Please check file naming or prefix.")

# ---------------- 5. Serial / Parallel Processing ----------------
def worker(tp):
    seg_p, bdg_p, out_p = tp
    print(f"[INFO] Mapping: {bdg_p.name}  ←  {seg_p.name}")
    map_seg.main(seg_p, bdg_p, out_p)

if args.threads > 1 and len(tasks) > 1:
    from multiprocessing import Pool
    with Pool(processes=args.threads) as pool:
        pool.map(worker, tasks)
else:
    for t in tasks:
        worker(t)

print(f"\n🎉 Successfully completed segmean mapping for {len(tasks)} bedGraph files. Results saved in: {args.out_dir}")