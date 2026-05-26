#!/bin/bash
# ======================================================
# Description: Extract FASTA sequences for INcen regions
# ======================================================

set -euo pipefail

# Path Configuration
BED_DIR="/share/org/YZWL/yzbsl_chengch/03.confirm_sat"
FA_DIR="/share/org/YZWL/yzbsl_chengch/13OC"
OUT_DIR="/share/org/YZWL/yzbsl_chengch/lxr/HOR/new/INcen_out"
mkdir -p "$OUT_DIR"

# Initialize Manifest File
MANIFEST="$OUT_DIR/incen_manifest.tsv"
echo -e "assembly\tbed_file\tfasta_file\tout_fasta\tstatus" > "$MANIFEST"

# Verify bedtools installation
if ! command -v bedtools >/dev/null 2>&1; then
  echo "[ERROR] bedtools is not installed or not found in PATH." >&2
  exit 1
fi

# Define supported FASTA file extensions
FA_EXTS=("fa" "fasta" "fna" "fa.gz" "fasta.gz" "fna.gz")

# ------------------------------------------------------
# Function: find_fasta_for
# Description: Locates matching FASTA files using direct 
#              matching followed by recursive fuzzy search.
# ------------------------------------------------------
find_fasta_for() {
  local key="$1"
  local fa=""
  
  # Attempt direct match with extensions
  for ext in "${FA_EXTS[@]}"; do
    if [[ -f "$FA_DIR/${key}.${ext}" ]]; then
      echo "$FA_DIR/${key}.${ext}"
      return 0
    fi
  done

  # Perform recursive fuzzy matching if direct match fails
  fa=$(find "$FA_DIR" -type f \
        \( -iname "${key}.fa*" -o -iname "${key}.fna*" \) 2>/dev/null | head -n 1)
  if [[ -n "$fa" ]]; then
    echo "$fa"
    return 0
  fi

  return 1
}

# ------------------------------------------------------
# Main Execution Loop
# ------------------------------------------------------
shopt -s nullglob
for bed in "$BED_DIR"/*.sat.sort.bed; do
  name=$(basename "$bed")
  key="${name%.sat.sort.bed}"
  out_fa="$OUT_DIR/${key}.INcen.fasta"

  echo ">>> [${key}] Processing..."

  # Locate the corresponding FASTA file
  fasta=$(find_fasta_for "$key" || true)
  if [[ -z "$fasta" ]]; then
    echo "[WARN] FASTA file not found for: $key" >&2
    echo -e "${key}\t${bed}\t\t${out_fa}\tMISSING_FASTA" >> "$MANIFEST"
    continue
  fi

  # Execute bedtools getfasta
  if bedtools getfasta -fi "$fasta" -bed "$bed" -fo "$out_fa" -name 2>/dev/null; then
    if [[ -s "$out_fa" ]]; then
      echo -e "${key}\t${bed}\t${fasta}\t${out_fa}\tOK" >> "$MANIFEST"
    else
      echo "[WARN] Output is empty: $out_fa (Possible chromosome name mismatch between BED and FASTA)" >&2
      echo -e "${key}\t${bed}\t${fasta}\t${out_fa}\tEMPTY" >> "$MANIFEST"
    fi
  else
    echo "[ERR] bedtools extraction failed for: $key" >&2
    echo -e "${key}\t${bed}\t${fasta}\t${out_fa}\tFAILED" >> "$MANIFEST"
  fi
done

echo
echo "[OK] All tasks completed. Manifest file generated at: $MANIFEST"