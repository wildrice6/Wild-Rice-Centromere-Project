#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HOR detection & scoring FROM precomputed pairwise identities.

Usage examples
--------------
# Single file
python hor_from_pairs.py \
  --fa /path/AA_Ogla_hap1.Chr01.w100.Top1.fa \
  --pairs /path/AA_Ogla_hap1.Chr01.w100.pairwise.tsv \
  --outdir results_pairs \
  --id_threshold 0.96 --scale auto \
  --min_block 3 --max_block 10 --max_shift 0

# Batch: many FASTAs in dir A, pairwise TSVs in dir B with matched basename + suffix
python hor_from_pairs.py \
  --fa_glob "/share/.../06.repeat_top1_fa/*.fa" \
  --pairs_dir "/share/.../pairwise_from_moddot/" \
  --pairs_suffix ".pairwise.tsv" \
  --outdir /share/.../hor_from_pairs_out \
  --id_threshold 0.96 --scale auto \
  --min_block 3 --max_block 10 --max_shift 0
"""

import argparse
import csv
import glob
import os
import sys
from typing import Dict, List, Tuple, Optional

# ---------- FASTA ----------
try:
    from Bio import SeqIO
except Exception:
    SeqIO = None


def read_fa_order(fa_path: str) -> List[str]:
    """Return labels in the order they appear in FASTA."""
    if SeqIO is None:
        raise RuntimeError("Biopython is required. Install: pip install biopython")
    labels = [str(rec.id) for rec in SeqIO.parse(fa_path, "fasta")]
    if not labels:
        raise ValueError(f"Empty FASTA: {fa_path}")
    return labels


# ---------- PAIRWISE TSV/CSV ----------
def sniff_delim(path: str) -> str:
    """Very simple delimiter sniffer: automatically determine if the file is tab or comma delimited."""
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        head = f.readline()
    return "\t" if "\t" in head else ","


def pairs_identity_max(path: str, col_iden: str, delim: Optional[str]) -> float:
    """First pass: get max identity value to determine scale (percentage vs fraction)."""
    if delim is None:
        delim = sniff_delim(path)
    mx = -1.0
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        rdr = csv.DictReader(f, delimiter=delim)
        for r in rdr:
            try:
                val = float(r[col_iden])
                if val > mx:
                    mx = val
            except Exception:
                continue
    if mx < 0:
        raise ValueError(f"No numeric identity found in: {path}")
    return mx


def build_similarity_from_pairs(
    labels: List[str],
    pairs_path: str,
    col_a: str,
    col_b: str,
    col_iden: str,
    id_threshold: float,
    scale: str = "auto",
    max_shift: int = 0,
    delim: Optional[str] = None,
) -> List[List[bool]]:
    """
    Build boolean similarity matrix S (N x N) using precomputed identities.
    labels: array order (from FASTA)
    id_threshold: threshold in fraction (e.g., 0.96). If scale=='percent', we compare to 96.
    scale: 'auto'|'percent'|'fraction'
    max_shift: 0 means no limit; else only |i-j|<=max_shift pairs are considered.
    """
    n = len(labels)
    idx = {lab: i for i, lab in enumerate(labels)}
    S = [[False] * n for _ in range(n)]
    for i in range(n):
        S[i][i] = True

    if delim is None:
        delim = sniff_delim(pairs_path)

    # Decide scale
    if scale == "auto":
        mx = pairs_identity_max(pairs_path, col_iden, delim)
        use_fraction = mx <= 1.5  # <=1.5 -> very likely fraction [0,1]
    elif scale == "fraction":
        use_fraction = True
    else:
        use_fraction = False  # percent

    thr = id_threshold if use_fraction else id_threshold * 100.0

    # Second pass: fill matrix
    miss_a = miss_b = 0
    total = 0
    with open(pairs_path, "r", encoding="utf-8", errors="ignore") as f:
        rdr = csv.DictReader(f, delimiter=delim)
        for r in rdr:
            total += 1
            try:
                a = r[col_a].strip()
                b = r[col_b].strip()
                val = float(r[col_iden])
            except Exception:
                continue
            ia = idx.get(a)
            if ia is None:
                miss_a += 1
                continue
            ib = idx.get(b)
            if ib is None:
                miss_b += 1
                continue
            if ia == ib:
                continue
            if max_shift > 0 and abs(ib - ia) > max_shift:
                continue
            if val >= thr:
                S[ia][ib] = True
                S[ib][ia] = True

    if (miss_a or miss_b) and total:
        print(f"[WARN] {os.path.basename(pairs_path)}: skipped {miss_a}+{miss_b} rows due to unmatched labels.", file=sys.stderr)
    return S


# ---------- HOR detection & scoring ----------
def find_hors(arr: List[str], S: List[List[bool]], min_block: int, max_block: int):
    """
    Slide offsets d = 1..N-1; a contiguous diagonal run of True of length Lrun
    yields HORs with block_size L where min_block <= L <= min(Lrun, max_block).
    """
    N = len(arr)
    hors = []
    for d in range(1, N):
        i = 0
        end = N - d
        while i < end:
            if S[i][i + d]:
                start = i
                Lrun = 1
                k = i + 1
                while k < end and S[k][k + d]:
                    Lrun += 1
                    k += 1
                maxL = min(Lrun, max_block)
                for L in range(min_block, maxL + 1):
                    for off in range(0, Lrun - L + 1):
                        a0 = start + off
                        b0 = a0 + d
                        Apos = list(range(a0 + 1, a0 + L + 1))
                        Bpos = list(range(b0 + 1, b0 + L + 1))
                        hors.append(
                            {
                                "block_size": L,
                                "shift": d,
                                "a_start": a0,
                                "b_start": b0,
                                "A_pos": Apos,
                                "B_pos": Bpos,
                                "A_labels": [arr[x - 1] for x in Apos],
                                "B_labels": [arr[x - 1] for x in Bpos],
                            }
                        )
                i = k
            else:
                i += 1
    return hors


def compute_scores(arr: List[str], hors: List[Dict]):
    """
    For each position i:
      forms[i]   = how many HORs the repeat at position i forms
      partners[i]= distinct labels it pairs with in the other block (symmetrical)
      score%     = |partners[i]| / N * 100
    """
    N = len(arr)
    forms = [0] * N
    partners = [set() for _ in range(N)]
    for h in hors:
        L = h["block_size"]
        a0 = h["a_start"]
        b0 = h["b_start"]
        for t in range(L):
            i = a0 + t
            j = b0 + t
            
            # --- [Core Modification] ---
            # Acknowledge that pairings are bidirectional
            
            # Unit i participates in an HOR
            forms[i] += 1
            # The partner of unit i is j
            partners[i].add(arr[j])
            
            # Unit j also participates in this HOR
            forms[j] += 1
            # The partner of unit j is i
            partners[j].add(arr[i])
            # --- [End Modification] ---
            
    rows = []
    for i, lab in enumerate(arr):
        uniq = len(partners[i])
        score_pct = 100.0 * uniq / N if N > 0 else 0.0
        rows.append(
            {
                "pos": i + 1,
                "label": lab,
                "HORs_the_repeat_forms": forms[i],
                "Unique_partners_in_2nd_block": uniq,
                "HOR_score_percent": f"{score_pct:.6f}",
            }
        )
    return rows

def write_pairs(path: str, hors: List[Dict]):
    with open(path, "w", encoding="utf-8") as f:
        f.write("HOR_ID\tblock_size\tshift\ta_start\tb_start\tA_pos\tB_pos\tA_labels\tB_labels\n")
        for k, h in enumerate(hors, 1):
            f.write(
                f"{k}\t{h['block_size']}\t{h['shift']}\t{h['a_start']+1}\t{h['b_start']+1}\t"
                f"{';'.join(map(str, h['A_pos']))}\t{';'.join(map(str, h['B_pos']))}\t"
                f"{';'.join(h['A_labels'])}\t{';'.join(h['B_labels'])}\n"
            )


def write_scores(path: str, rows: List[Dict]):
    with open(path, "w", encoding="utf-8") as f:
        f.write("pos\tlabel\tHORs_the_repeat_forms\tUnique_partners_in_2nd_block\tHOR_score_percent\n")
        for r in rows:
            f.write(
                f"{r['pos']}\t{r['label']}\t{r['HORs_the_repeat_forms']}\t"
                f"{r['Unique_partners_in_2nd_block']}\t{r['HOR_score_percent']}\n"
            )


# ---------- Runner ----------
def run_one(
    fa_path: str,
    pairs_path: str,
    outdir: str,
    id_threshold: float,
    min_block: int,
    max_block: int,
    max_shift: int,
    col_a: str,
    col_b: str,
    col_id: str,
    scale: str,
):
    base = os.path.splitext(os.path.basename(fa_path))[0]
    os.makedirs(outdir, exist_ok=True)
    out_pairs = os.path.join(outdir, base + ".hor_pairs.tsv")
    out_score = os.path.join(outdir, base + ".hor_score.tsv")

    arr = read_fa_order(fa_path)
    S = build_similarity_from_pairs(
        labels=arr,
        pairs_path=pairs_path,
        col_a=col_a,
        col_b=col_b,
        col_iden=col_id,
        id_threshold=id_threshold,
        scale=scale,
        max_shift=max_shift,
    )
    hors = find_hors(arr, S, min_block, max_block)
    rows = compute_scores(arr, hors)

    write_pairs(out_pairs, hors)
    write_scores(out_score, rows)
    print(f"[OK] {base}: N={len(arr)}, HORs={len(hors)} → {out_pairs} ; {out_score}")


def main():
    ap = argparse.ArgumentParser("HOR detection from precomputed pairwise identities")
    # single file mode
    ap.add_argument("--fa", help="FASTA file (array order)")
    ap.add_argument("--pairs", help="pairwise table for the FASTA")
    # batch mode
    ap.add_argument("--fa_glob", help="glob pattern for multiple FASTAs, e.g. '/dir/*.fa'")
    ap.add_argument(
        "--pairs_dir",
        help="directory holding pairwise TSVs (used with --fa_glob). "
             "Pairs path = pairs_dir / (basename(fa).replace('.fa', pairs_suffix))",
    )
    ap.add_argument(
        "--pairs_suffix",
        default=".pairwise.tsv",
        help="suffix to replace '.fa' when searching the pairwise file (default: .pairwise.tsv)",
    )
    ap.add_argument("--outdir", default="results_pairs")

    # pairwise columns & scale
    ap.add_argument("--col_a", default="win_id_a")
    ap.add_argument("--col_b", default="win_id_b")
    ap.add_argument("--identity_col", default="identity")
    ap.add_argument("--scale", choices=["auto", "percent", "fraction"], default="auto")

    # HOR params
    ap.add_argument("--id_threshold", type=float, default=0.96)
    ap.add_argument("--min_block", type=int, default=3)
    ap.add_argument("--max_block", type=int, default=10)
    ap.add_argument(
        "--max_shift",
        type=int,
        default=0,
        help="0=no limit; else only consider pairs with |i-j|<=max_shift",
    )

    args = ap.parse_args()

    if args.fa_glob:
        fas = sorted(glob.glob(args.fa_glob))
        if not fas:
            print(f"[ERR] No FASTA matched: {args.fa_glob}")
            sys.exit(1)
        for fa in fas:
            base = os.path.basename(fa).replace(".fa", "")
            if args.pairs_dir:
                pairs = os.path.join(args.pairs_dir, base + args.pairs_suffix)
            else:
                pairs = os.path.splitext(fa)[0] + args.pairs_suffix
            if not os.path.exists(pairs):
                print(f"[WARN] pairs not found for {fa}: {pairs} (skip)")
                continue
            run_one(
                fa_path=fa,
                pairs_path=pairs,
                outdir=args.outdir,
                id_threshold=args.id_threshold,
                min_block=args.min_block,
                max_block=args.max_block,
                max_shift=args.max_shift,
                col_a=args.col_a,
                col_b=args.col_b,
                col_id=args.identity_col,
                scale=args.scale,
            )
    else:
        if not args.fa or not args.pairs:
            print("[ERR] Provide --fa and --pairs, or use --fa_glob for batch mode.")
            sys.exit(1)
        run_one(
            fa_path=args.fa,
            pairs_path=args.pairs,
            outdir=args.outdir,
            id_threshold=args.id_threshold,
            min_block=args.min_block,
            max_block=args.max_block,
            max_shift=args.max_shift,
            col_a=args.col_a,
            col_b=args.col_b,
            col_id=args.identity_col,
            scale=args.scale,
        )


if __name__ == "__main__":
    main()