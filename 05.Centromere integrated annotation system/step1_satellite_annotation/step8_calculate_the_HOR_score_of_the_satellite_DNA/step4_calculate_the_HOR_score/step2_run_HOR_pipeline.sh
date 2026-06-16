#!/bin/bash -l
#CSUB -J HOR_Score_ALL
#CSUB -q c01
#CSUB -o /share/org/YZWL/yzbsl_chengch/lxr/HOR/job.o
#CSUB -e /share/org/YZWL/yzbsl_chengch/lxr/HOR/job.e
#CSUB -n 8
#CSUB -R "span[hosts=1]"
#CSUB -cwd /share/org/YZWL/yzbsl_chengch/lxr/HOR/new

set -euo pipefail

echo "=========================================================="
echo "JobID: ${LSB_JOBID:-NA}"
echo "Host : $(hostname)"
echo "CWD  : $(pwd)"
date
echo "=========================================================="

# Path Configuration
PYTHON_BIN="$HOME/.conda/envs/HOR_env/bin/python"
SCRIPT_HOR="/share/org/YZWL/yzbsl_chengch/lxr/HOR/CC126-155/Calculate_the_HOR_score.py"
STRICT_ROOT="/share/org/YZWL/yzbsl_chengch/lxr/HOR/new/INcen_labeled"
OUTROOT="/share/org/YZWL/yzbsl_chengch/lxr/HOR/new/HOR_results_all"
LOG_FAIL="${OUTROOT}/HOR_failed.log"

mkdir -p "$OUTROOT"
: > "$LOG_FAIL"  # Clear previous error logs

# Dependency Verification
[[ -x "$PYTHON_BIN" ]] || { echo "[FATAL] Python interpreter not found: $PYTHON_BIN"; exit 127; }
[[ -s "$SCRIPT_HOR"  ]] || { echo "[FATAL] Script missing: $SCRIPT_HOR"; exit 2; }

echo ">>> Utilizing Python environment: $PYTHON_BIN"
"$PYTHON_BIN" - <<'PY'
import pandas, openpyxl, sys
print(f"[Dependencies Verified] pandas version: {pandas.__version__}")
PY

# ------------------------------------------------------
# Function: fix_strict_if_needed
# Description: Resolves missing haplotype information in 
#              sequence identifiers within metadata files.
# ------------------------------------------------------
fix_strict_if_needed() {
  local strict_xlsx="$1"
  local assembly="$2"
  local cen="$3"
  local out_fixed="${strict_xlsx%.xlsx}.fixed.xlsx"

  "$PYTHON_BIN" - <<PY || { echo "$strict_xlsx"; exit 0; }
import re, sys
import pandas as pd
from pathlib import Path
src = Path("${strict_xlsx}")
dst = Path("${out_fixed}")
assembly = "${assembly}"
cen = "${cen}"
df = pd.read_excel(src)

# Validate mandatory columns
need_cols = {"id","code","length","type","chr"}
if not need_cols.issubset(df.columns):
    print(src)
    sys.exit(0)

# Determine haplotype designation
has_hap = bool(re.search(r"_hap\\d+$", assembly))
assembly_with_hap = assembly if has_hap else (assembly + "_hap1")

# Identifier reconstruction logic
pat = re.compile(r"^(CEN\\d+)_([A-Za-z0-9_]+)_(Chr\\d{2})_INcen_(\\d+)$")
def rebuild(row):
    s = str(row["id"])
    m = pat.match(s)
    if m:
        cen_t, asm, chr_s, idx = m.groups()
        if not re.search(r"_hap\\d+$", asm):
            asm = assembly_with_hap
        return f"{cen_t}_{asm}_{chr_s}_INcen_{idx}"
    
    m2 = re.search(r"_INcen_(\\d+)$", s)
    idx = m2.group(1) if m2 else None
    chr_s = row.get("chr")
    cen_t = str(row.get("type","")).upper() or cen
    if isinstance(chr_s, str) and re.match(r"Chr\\d{2}$", chr_s) and idx is not None:
        return f"{cen_t}_{assembly_with_hap}_{chr_s}_INcen_{idx}"
    return s

df["id"] = df.apply(rebuild, axis=1)
df.to_excel(dst, index=False)
print(dst)
PY
}

# ------------------------------------------------------
# Function: run_one
# Description: Processes a single assembly/CEN dataset.
# ------------------------------------------------------
run_one() {
  local strict_xlsx="$1"
  local cen="$(basename "$(dirname "$strict_xlsx")")"
  local assembly="$(basename "$(dirname "$(dirname "$strict_xlsx")")")"
  local hap="hapX"
  [[ "$assembly" =~ _hap([0-9]+)$ ]] && hap="hap${BASH_REMATCH[1]}"
  local outdir="${OUTROOT}/${assembly}/${cen}/strict_096_b3"

  echo ">>> Processing Dataset: ${assembly} / ${cen}"
  mkdir -p "$outdir"

  # Standardize identifiers missing haplotype designations
  fixed_path="$(fix_strict_if_needed "$strict_xlsx" "$assembly" "$cen" 2>/dev/null || echo "$strict_xlsx")"

  # Verify input file validity
  [[ -s "$fixed_path" ]] || { echo "[SKIP] File empty or missing: $fixed_path" | tee -a "$LOG_FAIL"; return; }

  # Execute Higher-Order Repeat (HOR) computation
  if ! "$PYTHON_BIN" "$SCRIPT_HOR" \
    --label_map "$fixed_path" \
    --id_col id \
    --code_col code \
    --outdir "$outdir" \
    --id_threshold 0.96 \
    --min_block 3 \
    --max_block 50 \
    --label_mode \
    --filter_cen "$cen" \
    --print_ids; then
      echo "[FAIL] Error encountered for $assembly / $cen using $fixed_path" | tee -a "$LOG_FAIL"
      return
  fi
}

# ------------------------------------------------------
# Main Execution: Traverse all standardized metadata files
# ------------------------------------------------------
mapfile -t STRICT_LIST < <(find "$STRICT_ROOT" -type f -name "*_STRICT.xlsx" | sort)
echo ">>> Total standardized metadata files to process: ${#STRICT_LIST[@]}"

for x in "${STRICT_LIST[@]}"; do
  echo "----------------------------------------------------------"
  run_one "$x" || true   # Continue execution even if a single task fails
done

echo "----------------------------------------------------------"
echo ">>> Pipeline execution finalized."
echo ">>> Output Directory: $OUTROOT"
echo ">>> Error Log: $LOG_FAIL"
echo "----------------------------------------------------------"