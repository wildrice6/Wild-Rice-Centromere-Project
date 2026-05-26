#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Batch convert 'label_map.xlsx' files located in /.../INcen_labeled/{assembly}/{CEN}/ 
into a strictly standardized format (columns: id, code, length, type, chr).

ID Convention: {type}_{assembly}_{chr}_INcen_{idx}

Key Improvements:
- Automatic Index Generation: If the 'id' in label_map does not contain "idx=...", 
  indices (idx=1..N) are automatically generated per chromosome based on coordinates 
  or original row order.
- Vectorized Processing: IDs are generated using vectorized operations to ensure 
  consistency and prevent length mismatch errors.
- Coordinate Compatibility: Falls back to original row sequence if genomic 
  coordinates are missing.
- Reference Synchronization: Optional support for a reference table 
  ({assembly}_{CEN}_STRICT.ref.xlsx) to override id/code based on 
  (chr, idx, length, type) matching.
"""

import re
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# ========= Configuration =========
ROOT = Path("/share/org/YZWL/yzbsl_chengch/lxr/HOR/new/INcen_labeled")  # Input root directory
LABEL_MAP_CANDIDATES = ("*_label_map.xlsx", "label_map.xlsx")           # Possible label_map filenames
OUT_NAME_TEMPLATE = "{assembly}_{cen}_STRICT.xlsx"                      # Output filename template
REF_NAME_TEMPLATE = "{assembly}_{cen}_STRICT.ref.xlsx"                  # Optional reference table template

def extract_chr(s: str) -> str | None:
    """Extracts chromosome name and standardizes it to 'ChrXX' format."""
    m = re.search(r"(Chr\d{1,2}|chr\d{1,2})", str(s))
    if not m:
        return None
    num = re.search(r"\d{1,2}", m.group(1))
    return f"Chr{int(num.group(0)):02d}" if num else m.group(1)

def extract_idx_from_label_id(s: str) -> float | None:
    """Extracts index from strings containing 'idx=NNN'."""
    m = re.search(r"idx=(\d+)", str(s))
    return float(m.group(1)) if m else np.nan

def extract_start_pos(s: str) -> float | None:
    """Extracts the start coordinate from strings like 'Chr01:12345-67890'."""
    m = re.search(r":(\d+)-(\d+)", str(s))
    return float(m.group(1)) if m else np.nan

def parse_idx_from_strict_id(s: str) -> float | None:
    """Extracts index from strict IDs ending in '_INcen_NNN'."""
    m = re.search(r"_INcen_(\d+)$", str(s))
    return float(m.group(1)) if m else np.nan

def find_label_map_file(cen_dir: Path) -> Path | None:
    """Locates the label_map.xlsx file within a specific CEN directory."""
    for pat in LABEL_MAP_CANDIDATES:
        hits = sorted(cen_dir.glob(pat))
        if hits:
            return hits[0]
    return None

def main():
    if not ROOT.exists():
        print(f"[ERROR] Root directory does not exist: {ROOT}", file=sys.stderr)
        sys.exit(1)

    total_dirs = 0
    made_files = 0
    skipped = 0

    # Iterate through assembly directories
    for assembly_dir in sorted(p for p in ROOT.iterdir() if p.is_dir()):
        assembly = assembly_dir.name  # e.g., CC_Ooff_hap2
        
        # Iterate through CEN type directories
        for cen_dir in sorted(p for p in assembly_dir.iterdir() if p.is_dir()):
            cen = cen_dir.name        # e.g., CEN126 / CEN155
            total_dirs += 1

            lm = find_label_map_file(cen_dir)
            if lm is None:
                print(f"[WARN] No label_map.xlsx found in: {cen_dir}")
                skipped += 1
                continue

            try:
                df = pd.read_excel(lm, sheet_name=0)
            except Exception as e:
                print(f"[ERR] Failed to read file: {lm} ({e})")
                skipped += 1
                continue

            # Verify required columns
            required_cols = {"id", "code", "length", "type"}
            if not required_cols.issubset(df.columns):
                print(f"[ERR] Missing required columns {required_cols} in: {lm}")
                skipped += 1
                continue

            # —— Parse Chromosome (chr)
            if "chr" in df.columns and df["chr"].notna().any():
                # If chr column exists, standardize to ChrXX
                chr_guess = df["chr"].astype(str).map(extract_chr)
            else:
                chr_guess = df["id"].astype(str).map(extract_chr)

            # —— Parse Index (idx)
            # Preference: 1. idx=NNN in id string | 2. Automatic numbering
            idx_guess = df["id"].astype(str).map(extract_idx_from_label_id)

            # —— Handle missing indices: Auto-numbering based on chr/coordinates/row order
            start_guess = df["id"].astype(str).map(extract_start_pos)
            order_col = np.arange(len(df))  # Preserve original row order

            tmp = pd.DataFrame({
                "_chr": chr_guess,
                "_idx": idx_guess,
                "_start": start_guess,
                "_order": order_col
            })

            miss_chr = tmp["_chr"].isna().sum()
            miss_idx = tmp["_idx"].isna().sum()
            if miss_chr > 0 or miss_idx > 0:
                print(f"[WARN] {assembly}/{cen} parsing issues: missing {miss_chr} chr, {miss_idx} idx. Applying auto-numbering.")

            # Calculate natural order per chromosome: prioritize start pos, then row order
            sorter = tmp[["_chr", "_start", "_order"]].copy()
            sorter["_rank"] = sorter.groupby("_chr", dropna=False)["_start"] \
                                   .rank(method="first", na_option="bottom")
            
            # Sort indices for group processing
            sort_index = sorter.sort_values(by=["_chr", "_rank", "_order"], kind="mergesort").index

            # Calculate cumulative count (1..N) within each chromosome group
            idx_auto = pd.Series(index=tmp.index, dtype="float64")
            sorted_tmp = tmp.loc[sort_index, :]
            auto_seq = sorted_tmp.groupby("_chr", dropna=False).cumcount() + 1
            idx_auto.loc[sort_index] = auto_seq.values
            
            # Fill missing indices with generated values
            idx_final = tmp["_idx"].copy()
            idx_final = idx_final.where(idx_final.notna(), idx_auto)

            # Drop records that still lack chr or idx information
            keep = idx_final.notna() & tmp["_chr"].notna()
            if keep.sum() < len(df):
                drop_n = int((~keep).sum())
                print(f"[WARN] {assembly}/{cen}: {drop_n} rows discarded due to missing chr/idx.")
            
            df2 = df.loc[keep, ["id", "code", "length", "type"]].copy()
            chr_final = tmp.loc[keep, "_chr"].astype(str)
            idx_final = idx_final.loc[keep].astype("int64")

            # —— Construct Strict ID (Vectorized)
            # Format: {type}_{assembly}_{chr}_INcen_{idx}
            id_strict = (
                df2["type"].astype(str).str.upper() + "_" +
                assembly + "_" +
                chr_final + "_INcen_" +
                idx_final.astype(str)
            )

            strict_df = pd.DataFrame({
                "id": id_strict.values,
                "code": df2["code"].astype(str).values,
                "length": df2["length"].values,
                "type": df2["type"].values,
                "chr": chr_final.values
            })

            # —— Optional: Reference Synchronization
            ref_path = cen_dir / REF_NAME_TEMPLATE.format(assembly=assembly, cen=cen)
            if ref_path.exists():
                try:
                    ref = pd.read_excel(ref_path, sheet_name=0)
                    ref = ref.copy()
                    ref["idx"] = ref["id"].astype(str).map(parse_idx_from_strict_id)
                    
                    strict_df = strict_df.copy()
                    strict_df["idx"] = strict_df["id"].astype(str).map(parse_idx_from_strict_id)

                    merged = pd.merge(
                        ref[["chr","idx","length","type","id","code"]].rename(columns={"id":"id_ref","code":"code_ref"}),
                        strict_df[["chr","idx","length","type","id","code"]],
                        on=["chr","idx","length","type"],
                        how="right",
                        validate="one_to_one"
                    )
                    merged["id"]   = merged["id_ref"].combine_first(merged["id"])
                    merged["code"] = merged["code_ref"].combine_first(merged["code"])
                    strict_df = merged[["id","code","length","type","chr"]].copy()
                    print(f"[INFO] {assembly}/{cen}: id/code overridden by reference table: {ref_path.name}")
                except Exception as e:
                    print(f"[WARN] Override failed (ignoring reference): {ref_path} ({e})")

            # —— Output Generation
            out_path = cen_dir / OUT_NAME_TEMPLATE.format(assembly=assembly, cen=cen)
            try:
                strict_df.to_excel(out_path, index=False)
                made_files += 1
                print(f"[OK] Output generated: {out_path} (Rows: {len(strict_df)})")
            except Exception as e:
                print(f"[ERR] Failed to write output: {out_path} ({e})")
                skipped += 1

    print(f"\n[SUMMARY] Directories processed: {total_dirs} | Successful: {made_files} | Skipped/Failed: {skipped}")

if __name__ == "__main__":
    main()