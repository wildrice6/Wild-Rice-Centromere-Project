#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description:
Batch process extracted INcen FASTA files (supporting .fa/.fasta/.fna and .gz) from a directory.
Sequences are categorized by (assembly, CEN type) and processed as follows:
  1) Generate a fully labeled FASTA: >{CODE}|{OriginalID}
  2) Generate a deduplicated FASTA containing only unique sequences.
  3) Generate a mapping table (XLSX): OriginalID → CODE, length, type, assembly, source file.

Output Directory Structure:
/share/.../INcen_labeled/
└─ {assembly}/
   └─ {CEN_type}/
      ├─ {assembly}_{CEN_type}_labeled.full.fa
      ├─ {assembly}_{CEN_type}_labeled.unique.fa
      └─ {assembly}_{CEN_type}_label_map.xlsx

Notes:
- CEN type identification is not limited to CEN126/155; it matches patterns like r"CEN\\d+" (case-insensitive).
- The assembly name is primarily inferred from the filename (e.g., CC_Ooff_hap1) and secondarily from the sequence header.
- Deduplication and labeling (A1, B1, C1...) are performed independently for each (assembly, CEN) group.
"""

import os
import re
import sys
import gzip
import hashlib
from pathlib import Path
from typing import Iterator, Tuple, Dict, List
import pandas as pd

# ========= Configuration =========
INPUT_DIR  = "/share/org/YZWL/yzbsl_chengch/lxr/HOR/new/INcen_out"
OUTPUT_DIR = "/share/org/YZWL/yzbsl_chengch/lxr/HOR/new/INcen_labeled"

ACCEPT_EXTS = {".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"}
WRAP_WIDTH  = 80
USE_UPPERCASE_CODES = True  # Use uppercase prefixes: A1, B1, ...

# ========= Utility Functions =========
def open_maybe_gz(path: str):
    """Opens a file, handling Gzip compression if the suffix is .gz."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def fasta_reader(path: str) -> Iterator[Tuple[str, str]]:
    """Simple FASTA reader that yields (header_without_>, sequence)."""
    header = None
    seq_chunks = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)

def wrap_sequence(seq: str, width: int = 80) -> str:
    """Wraps FASTA sequences to a specified character width."""
    if width and width > 0:
        return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))
    else:
        return seq

def letter_labels(uppercase: bool = False) -> Iterator[str]:
    """Generates alphabetic labels (A, B, C... Z, AA, AB...). Used for sequence coding."""
    import string
    letters = string.ascii_uppercase if uppercase else string.ascii_lowercase
    n = 1
    while True:
        total = len(letters) ** n
        for i in range(total):
            s = ""
            x = i
            for _ in range(n):
                s = letters[x % len(letters)] + s
                x //= len(letters)
            yield s
        n += 1

def sha256_of_seq(seq: str) -> str:
    """Generates a SHA256 hash for a given sequence to facilitate deduplication."""
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()

def detect_cen_type(text: str) -> str:
    """Identifies the CEN type from text (e.g., CEN126, CEN155). Returns uppercase string or empty."""
    m = re.search(r"CEN\s*([0-9]+)", text, flags=re.IGNORECASE)
    if m:
        return f"CEN{m.group(1)}".upper()
    return ""

def detect_assembly_from_filename(name: str) -> str:
    """
    Infers the assembly name from the filename.
    Example Logic:
      CC_Ooff_hap1.INcen.CEN126.fasta  -> CC_Ooff_hap1
      AA_Osat_hap2.INcen.fa           -> AA_Osat_hap2
    """
    stem = name
    for sep in ["_CEN", ".CEN", "_cen", ".cen"]:
        if sep in stem:
            stem = stem.split(sep)[0]
            break
    for sep in ["_INcen", ".INcen", "_incen", ".incen"]:
        if sep in stem:
            stem = stem.split(sep)[0]
            break
    # Remove extensions
    stem = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", stem, flags=re.IGNORECASE)
    # Match standard format: Xxx_Yyy_hapN
    m = re.match(r"([A-Za-z0-9]+_[A-Za-z0-9]+_hap\d+)", stem)
    if m:
        return m.group(1)
    # Search for any string ending in hapN
    m2 = re.search(r"(.*?_hap\d+)", stem)
    if m2:
        return m2.group(1)
    # Fallback: split by delimiters
    return re.split(r"[.\s-]", stem)[0]

def detect_assembly_from_header(header: str) -> str:
    """Attempts to extract assembly name from the sequence header."""
    m = re.search(r"([A-Za-z0-9]+_[A-Za-z0-9]+_hap\d+)", header)
    if m:
        return m.group(1)
    return ""

def ensure_dir(p: Path):
    """Ensures that the target directory exists."""
    p.mkdir(parents=True, exist_ok=True)

# ========= Main Pipeline =========
def main():
    in_dir = Path(INPUT_DIR)
    out_root = Path(OUTPUT_DIR)
    ensure_dir(out_root)

    # Collect input FASTA files
    inputs = []
    for p in in_dir.rglob("*"):
        if not p.is_file():
            continue
        suf = "".join(p.suffixes).lower()
        if any(suf.endswith(ext) for ext in ACCEPT_EXTS):
            inputs.append(p)
    
    if not inputs:
        print(f"[ERROR] No FASTA files found in {in_dir} (Supported: {', '.join(ACCEPT_EXTS)})", file=sys.stderr)
        sys.exit(1)

    # Partition data into (assembly, cen_type) buckets
    buckets: Dict[Tuple[str,str], List[Tuple[str,str,str]]] = {}
    # Track skipped sequences
    skipped = []

    # Read files and categorize
    for fpath in sorted(inputs):
        fname = fpath.name
        file_cen = detect_cen_type(fname)
        file_assembly = detect_assembly_from_filename(fname)

        for header, seq in fasta_reader(str(fpath)):
            cen = detect_cen_type(header) or file_cen
            asm = detect_assembly_from_header(header) or file_assembly

            if not cen:
                skipped.append((str(fpath), header, "NO_CEN_TYPE"))
                continue
            if not asm:
                skipped.append((str(fpath), header, "NO_ASSEMBLY"))
                continue

            key = (asm, cen)
            buckets.setdefault(key, []).append((header, seq, str(fpath)))

    # Process each bucket
    summary_rows = []

    for (assembly, cen_type), items in sorted(buckets.items()):
        out_dir = out_root / assembly / cen_type
        ensure_dir(out_dir)

        # Output filenames
        base = f"{assembly}_{cen_type}_labeled"
        full_fa = out_dir / f"{base}.full.fa"
        uniq_fa = out_dir / f"{base}.unique.fa"
        map_xlsx = out_dir / f"{assembly}_{cen_type}_label_map.xlsx"

        # Deduplication and Labeling
        seqhash2code: Dict[str,str] = {}
        code2seq: Dict[str,str] = {}
        labeler = letter_labels(uppercase=USE_UPPERCASE_CODES)

        map_rows = []  # columns: id, code, length, type, assembly, source_file
        total = 0

        with open(full_fa, "w") as wfull:
            for header, seq, src in items:
                total += 1
                h = sha256_of_seq(seq)
                if h not in seqhash2code:
                    code = next(labeler) + "1"
                    seqhash2code[h] = code
                    code2seq[code] = seq
                else:
                    code = seqhash2code[h]
                
                map_rows.append((header, code, len(seq), cen_type, assembly, src))
                wfull.write(f">{code}|{header}\n")
                wfull.write(wrap_sequence(seq, WRAP_WIDTH) + "\n")

        with open(uniq_fa, "w") as wuniq:
            for code, s in code2seq.items():
                wuniq.write(f">{code}\n")
                wuniq.write(wrap_sequence(s, WRAP_WIDTH) + "\n")

        # Save mapping table to Excel
        df = pd.DataFrame(map_rows, columns=["id","code","length","type","assembly","source_file"])
        with pd.ExcelWriter(map_xlsx, engine="openpyxl") as xw:
            df.to_excel(xw, index=False, sheet_name="map")

        unique_n = len(code2seq)
        print(f"[OK] {assembly} / {cen_type}: total={total}, unique={unique_n}")
        print(f"     - Labeled FASTA: {full_fa}")
        print(f"     - Unique FASTA:  {uniq_fa}")
        print(f"     - Map Table:     {map_xlsx}")

        summary_rows.append((assembly, cen_type, total, unique_n, str(out_dir)))

    # Generate summary manifest
    summary_tsv = Path(OUTPUT_DIR) / "summary_manifest.tsv"
    if summary_rows:
        sdf = pd.DataFrame(summary_rows, columns=["assembly","type","total","unique","out_dir"])
        sdf.sort_values(["assembly","type"], inplace=True)
        sdf.to_csv(summary_tsv, sep="\t", index=False)
        print(f"[SUMMARY] -> {summary_tsv}")

    # Log skipped sequences if any
    if skipped:
        skipped_tsv = Path(OUTPUT_DIR) / "skipped_sequences.tsv"
        kdf = pd.DataFrame(skipped, columns=["file","header","reason"])
        kdf.to_csv(skipped_tsv, sep="\t", index=False)
        print(f"[WARN] Some sequences were skipped. Details: {skipped_tsv}")

if __name__ == "__main__":
    main()