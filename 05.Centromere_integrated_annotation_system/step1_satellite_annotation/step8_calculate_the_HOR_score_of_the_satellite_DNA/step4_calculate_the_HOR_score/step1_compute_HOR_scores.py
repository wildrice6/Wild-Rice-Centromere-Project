#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Pipeline for computing Higher-Order Repeat (HOR) scores: 
Excel Metadata -> Per-Chromosome Array (+ FASTA) -> HOR Analysis.

Key Refinements:
- Resolved column name redundancy (specifically 'id') to prevent dictionary key errors 
  when handling Series objects.
- Standardized access to columns using args.id_col and args.code_col.
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Optional

import pandas as pd

try:
    from Bio import SeqIO
    from Bio.Align import PairwiseAligner
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except Exception:
    SeqIO = None
    PairwiseAligner = None


# ------------------------- ID Parser -------------------------
# Matches strings like: CEN126_Oryza_sativa_hap1_Chr01_INcen_001
ID_RE = re.compile(
    r"^(CEN\d+)_([A-Za-z]+_[A-Za-z]+)_(hap\d+)_Chr0*([0-9]+)_INcen_([0-9]+)$"
)

def parse_id_to_meta(idstr: str):
    """Extracts metadata from standardized sequence identifiers."""
    m = ID_RE.match(idstr)
    if not m:
        return None
    cen, species, hap, chrnum, index = m.groups()
    return {
        # Note: 'id' is excluded here to avoid namespace collisions with original columns
        "cen": cen,
        "species": species,
        "hap": hap,
        "chr": f"Chr{int(chrnum):02d}",
        "index": int(index),
    }


# ------------------------- Excel IO -------------------------
def read_label_map_xlsx(xlsx_path: str, id_col: str, code_col: str) -> pd.DataFrame:
    """Reads mapping metadata from Excel and validates required columns."""
    df = pd.read_excel(xlsx_path, engine="openpyxl")

    # Retain mandatory columns only
    for col in (id_col, code_col):
        if col not in df.columns:
            raise ValueError(f"Excel file must contain columns: '{id_col}' and '{code_col}'")
    df = df[[id_col, code_col]].dropna()

    # Standardize data types and trim whitespace
    df[id_col] = df[id_col].astype(str).str.strip()
    df[code_col] = df[code_col].astype(str).str.strip()

    # Parse metadata (excluding 'id' to avoid duplication)
    meta_series = df[id_col].apply(parse_id_to_meta)
    ok = meta_series.notna()
    if not ok.any():
        raise ValueError("No identifiers could be parsed. Please check the ID format and regex pattern.")
    meta_df = pd.DataFrame(list(meta_series[ok]))

    # Merge: Retain a single unique ID column along with extracted metadata
    df_ok = df.loc[ok].reset_index(drop=True)
    out = pd.concat([df_ok.reset_index(drop=True), meta_df.reset_index(drop=True)], axis=1)

    # Fallback: Ensure column uniqueness
    if out.columns.duplicated().any():
        out = out.loc[:, ~out.columns.duplicated()]
    
    return out  # Columns: id, code, cen, species, hap, chr, index


# ------------------------- FASTA IO -------------------------
def load_fasta_id2seq(path: str) -> Dict[str, str]:
    """Loads FASTA records into a dictionary mapping IDs to sequences."""
    if SeqIO is None:
        raise RuntimeError("Biopython is required. Install via: pip install biopython")
    id2seq = {}
    for rec in SeqIO.parse(path, "fasta"):
        id2seq[str(rec.id)] = str(rec.seq).upper()
    if not id2seq:
        raise ValueError(f"FASTA file is empty or invalid: {path}")
    return id2seq

def write_code_fasta(path: str, code2seq: Dict[str, str]):
    """Writes a FASTA file using codes as sequence identifiers."""
    if SeqIO is None:
        raise RuntimeError("Biopython is required.")
    records = [SeqRecord(Seq(seq), id=code, description="") for code, seq in code2seq.items()]
    SeqIO.write(records, path, "fasta")


# ------------------------- HOR core -------------------------
def letter_key(label: str) -> str:
    """Extracts the prefix key from a label (e.g., 'A' from 'A1')."""
    m = re.match(r"^([A-Za-z]+)", label)
    return m.group(1) if m else label

def make_aligner():
    """Initializes a PairwiseAligner with standardized alignment parameters."""
    if PairwiseAligner is None:
        return None
    al = PairwiseAligner()
    al.mode = "global"
    al.match_score = 1.0
    al.mismatch_score = -1.0
    al.open_gap_score = -2.0
    al.extend_gap_score = -1.0
    return al

def identity_align(seq1: str, seq2: str, al) -> float:
    """Calculates identity score between two sequences via global alignment."""
    aln = al.align(seq1, seq2)
    if not aln:
        return 0.0
    a = aln[0]
    L = max(len(seq1), len(seq2))
    if L == 0:
        return 0.0
    return max(0.0, min(1.0, a.score / L))

def build_similarity(arr: List[str],
                     seqs_by_label: Optional[Dict[str, str]],
                     id_thr: float,
                     max_shift: int,
                     label_mode: bool) -> List[List[bool]]:
    """Constructs a similarity matrix based on sequence identity or label matching."""
    N = len(arr)
    S = [[False]*N for _ in range(N)]
    for i in range(N):
        S[i][i] = True

    if label_mode or seqs_by_label is None:
        keys = [letter_key(x) for x in arr]
        for i in range(N):
            jmax = N if max_shift <= 0 else min(N, i + max_shift + 1)
            for j in range(i+1, jmax):
                S[i][j] = S[j][i] = (keys[i] == keys[j])
        return S

    al = make_aligner()
    if al is None:
        raise RuntimeError("Bio.Align.PairwiseAligner unavailable. Install biopython or use --label_mode.")

    for i in range(N):
        si = seqs_by_label.get(arr[i])
        if si is None:
            raise ValueError(f"Sequence missing in FASTA for label: {arr[i]}")
        jmax = N if max_shift <= 0 else min(N, i + max_shift + 1)
        for j in range(i+1, jmax):
            sj = seqs_by_label.get(arr[j])
            if sj is None:
                raise ValueError(f"Sequence missing in FASTA for label: {arr[j]}")
            S[i][j] = S[j][i] = (identity_align(si, sj, al) >= id_thr)
    return S

def find_hors(arr: List[str], S: List[List[bool]], min_block: int, max_block: int):
    """Identifies Higher-Order Repeat patterns within the sequence array."""
    N = len(arr)
    hors = []
    for d in range(1, N):
        i = 0
        end = N - d
        while i < end:
            if S[i][i+d]:
                start = i
                Lrun = 1
                k = i + 1
                while k < end and S[k][k+d]:
                    Lrun += 1
                    k += 1
                maxL = min(Lrun, max_block)
                for L in range(min_block, maxL+1):
                    for off in range(0, Lrun - L + 1):
                        a0 = start + off
                        b0 = a0 + d
                        Apos = list(range(a0+1, a0+L+1))
                        Bpos = list(range(b0+1, b0+L+1))
                        hors.append({
                            "block_size": L, "shift": d,
                            "a_start": a0, "b_start": b0,
                            "A_pos": Apos, "B_pos": Bpos,
                            "A_labels": [arr[x-1] for x in Apos],
                            "B_labels": [arr[x-1] for x in Bpos],
                        })
                i = k
            else:
                i += 1
    return hors

def compute_scores(arr: List[str], hors: List[Dict]):
    """Computes HOR occurrence counts and partner uniqueness scores for each position."""
    N = len(arr)
    forms = [0]*N
    partners = [set() for _ in range(N)]
    for h in hors:
        L = h["block_size"]; a0 = h["a_start"]; b0 = h["b_start"]
        for t in range(L):
            i = a0 + t
            j = b0 + t
            forms[i] += 1
            partners[i].add(arr[j])
    rows = []
    for i, lab in enumerate(arr):
        uniq = len(partners[i])
        score_pct = 100.0 * uniq / N if N > 0 else 0.0
        rows.append({
            "pos": i+1, "label": lab,
            "HORs_the_repeat_forms": forms[i],
            "Unique_partners_in_2nd_block": uniq,
            "HOR_score_percent": f"{score_pct:.6f}",
        })
    return rows

def write_pairs(path: str, hors: List[Dict]):
    """Outputs identified HOR pairs to a TSV file."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("HOR_ID\tblock_size\tshift\ta_start\tb_start\tA_pos\tB_pos\tA_labels\tB_labels\n")
        for k, h in enumerate(hors, 1):
            f.write(
                f"{k}\t{h['block_size']}\t{h['shift']}\t{h['a_start']+1}\t{h['b_start']+1}\t"
                f"{';'.join(map(str,h['A_pos']))}\t{';'.join(map(str,h['B_pos']))}\t"
                f"{';'.join(h['A_labels'])}\t{';'.join(h['B_labels'])}\n"
            )

def write_scores(path: str, rows: List[Dict]):
    """Outputs computed HOR scores to a TSV file."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("pos\tlabel\tHORs_the_repeat_forms\tUnique_partners_in_2nd_block\tHOR_score_percent\n")
        for r in rows:
            f.write(
                f"{r['pos']}\t{r['label']}\t{r['HORs_the_repeat_forms']}\t"
                f"{r['Unique_partners_in_2nd_block']}\t{r['HOR_score_percent']}\n"
            )


# ------------------------- Main -------------------------
def main():
    ap = argparse.ArgumentParser(description="Convert Excel and FASTA data into per-chromosome HOR arrays.")
    ap.add_argument("--label_map", required=True, help="Input Excel label map file.")
    ap.add_argument("--fasta_global", help="Global FASTA file (header must match ID).")
    ap.add_argument("--outdir", default="results", help="Directory for output results.")

    ap.add_argument("--id_col", default="id", help="Name of the ID column in Excel.")
    ap.add_argument("--code_col", default="code", help="Name of the code column in Excel.")

    ap.add_argument("--filter_species", default="", help="Comma-separated species filter.")
    ap.add_argument("--filter_haps", default="all", help="Comma-separated haplotype filter.")
    ap.add_argument("--filter_cen", default="", help="Comma-separated CEN type filter.")
    ap.add_argument("--filter_chrs", default="all", help="Comma-separated chromosome filter.")

    ap.add_argument("--id_threshold", type=float, default=0.96, help="Identity threshold for sequence matching.")
    ap.add_argument("--min_block", type=int, default=3, help="Minimum block size for HOR.")
    ap.add_argument("--max_block", type=int, default=10, help="Maximum block size for HOR.")
    ap.add_argument("--max_shift", type=int, default=400, help="Maximum shift distance.")

    ap.add_argument("--label_mode", action="store_true", help="Perform HOR analysis based on labels instead of sequences.")
    ap.add_argument("--print_ids", action="store_true", help="Print sample IDs matching the filters.")
    ap.add_argument("--dry_run", action="store_true", help="Perform setup without executing HOR computation.")

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = read_label_map_xlsx(args.label_map, args.id_col, args.code_col)

    def to_set(s): return set([t.strip() for t in s.split(",") if t.strip()]) if s else set()
    keep_species = to_set(args.filter_species)
    keep_cen     = to_set(args.filter_cen)
    keep_haps    = to_set("" if args.filter_haps == "all" else args.filter_haps)
    keep_chrs    = to_set("" if args.filter_chrs == "all" else args.filter_chrs)

    m = df.copy()
    if keep_species: m = m[m["species"].isin(keep_species)]
    if keep_cen:     m = m[m["cen"].isin(keep_cen)]
    if keep_haps:    m = m[m["hap"].isin(keep_haps)]
    if keep_chrs:    m = m[m["chr"].isin(keep_chrs)]
    
    if m.empty:
        print("[WARN] No IDs matched the criteria. Please verify your filter settings.")
        sys.exit(0)

    if args.print_ids:
        print("[INFO] Sample matched IDs (First 30):")
        for s in m[args.id_col].head(30):
            print("  ", s)

    id2seq = None
    if not args.label_mode:
        if not args.fasta_global:
            raise ValueError("Sequence mode requires --fasta_global. Alternatively, use --label_mode.")
        id2seq = load_fasta_id2seq(args.fasta_global)

    summary = []

    # Group processing by Species, Haplotype, CEN type, and Chromosome
    for (species, hap, cen, chrname), sub in m.groupby(["species","hap","cen","chr"], sort=True):
        sub = sub.sort_values("index", kind="stable")
        codes = sub[args.code_col].tolist()
        csv_rows = len(codes)

        prefix   = f"{species}_{hap}_{cen}_{chrname}"
        out_base = os.path.join(args.outdir, prefix)

        # Output array: One 'code' per record
        with open(out_base + ".array.txt", "w", encoding="utf-8") as g:
            g.write("\n".join(codes) + "\n")

        # Output FASTA: One representative sequence per unique code
        if not args.label_mode:
            uniq_codes = list(dict.fromkeys(codes))
            code2seq: Dict[str, str] = {}

            # Preferred representative: First ID encountered in the current group
            first_id_by_code = (
                sub.drop_duplicates(args.code_col)
                   .set_index(args.code_col)[args.id_col]
                   .to_dict()
            )

            for c in uniq_codes:
                seq = None
                # 1) Attempt group-specific lookup
                if c in first_id_by_code:
                    cand_id = first_id_by_code[c]
                    seq = id2seq.get(cand_id)

                # 2) Fallback: Species and CEN specific scope
                if seq is None:
                    sub2 = m[(m["species"] == species) & (m["cen"] == cen) & (m[args.code_col] == c)]
                    if not sub2.empty:
                        cand_id = sub2.iloc[0][args.id_col]
                        seq = id2seq.get(cand_id)

                # 3) Final Fallback: Search entire table
                if seq is None:
                    sub3 = df[df[args.code_col] == c]
                    if not sub3.empty:
                        cand_id = sub3.iloc[0][args.id_col]
                        seq = id2seq.get(cand_id)

                if seq is None:
                    print(f"[WARN] {prefix}: Representative sequence for code={c} not found in FASTA (skipping).")
                else:
                    code2seq[c] = seq

            if not code2seq:
                print(f"[SKIP] {prefix}: Per-Chromosome FASTA is empty.")
                continue

            write_code_fasta(out_base + ".fasta", code2seq)

        if args.dry_run:
            summary.append({
                "species": species, "hap": hap, "cen": cen, "chr": chrname,
                "csv_rows": csv_rows, "array_len": csv_rows,
                "hors_detected": "NA",
                "out_array": out_base + ".array.txt",
                "out_pairs": "NA",
                "out_score": "NA",
            })
            print(f"[DRY] {prefix}: Array length = {csv_rows}")
            continue

        # HOR Computation
        try:
            seqs_by_label = None if args.label_mode else code2seq
            S = build_similarity(codes, seqs_by_label, args.id_threshold, args.max_shift, args.label_mode)
            hors = find_hors(codes, S, args.min_block, args.max_block)
            rows = compute_scores(codes, hors)
        except Exception as e:
            print(f"[ERROR] Processing failed for {prefix}: {e}")
            continue

        write_pairs(out_base + ".hor_pairs.tsv", hors)
        write_scores(out_base + ".hor_score.tsv", rows)

        summary.append({
            "species": species, "hap": hap, "cen": cen, "chr": chrname,
            "csv_rows": csv_rows, "array_len": len(codes),
            "hors_detected": len(hors),
            "out_array": out_base + ".array.txt",
            "out_pairs": out_base + ".hor_pairs.tsv",
            "out_score": out_base + ".hor_score.tsv",
        })
        print(f"[OK] {prefix}: Array length = {csv_rows}, HORs detected = {len(hors)}")

    # Generate final summary report
    sum_path = os.path.join(args.outdir, "summary.tsv")
    with open(sum_path, "w", encoding="utf-8") as f:
        cols = ["species","hap","cen","chr","csv_rows","array_len","hors_detected","out_array","out_pairs","out_score"]
        f.write("\t".join(cols) + "\n")
        for r in summary:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"[DONE] Summary report generated at: {sum_path}")

if __name__ == "__main__":
    main()