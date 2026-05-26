#!/bin/bash
# Description: direct_merge_CEN155.sh - Direct extraction and merging of genomic intervals.

INPUT_DIR="/share/org/YZWL/yzwl_shahd/alltest/0823PAN_CEN/47.multi_omics/07.satellite"
OUTPUT_FILE="/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/all_CEN155_merged.tsv"

echo "Initiating extraction and consolidation of CEN155 genomic intervals..."
echo ""

# Create header for the output file
echo -e "Sample\tChromosome\tStart\tEnd" > "$OUTPUT_FILE"

total=0
for file in "$INPUT_DIR"/*.satellite.bed; do
    # Skip files prefixed with 'XX_'
    filename=$(basename "$file")
    if [[ "$filename" == XX_* ]]; then
        continue
    fi
    
    sample="${filename%.satellite.bed}"
    
    # Utilize awk to extract CEN155 regions and append directly to the output file
    awk -F'\t' -v s="$sample" '
    $4 == "CEN155" {
        print s "\t" $1 "\t" $2 "\t" $3;
        count++;
    }
    END {
        if (count > 0) {
            printf "  %s: %d CEN155 regions identified\n", s, count > "/dev/stderr";
        }
    }
    ' "$file" >> "$OUTPUT_FILE"
    
    # Calculate count for final summary
    count=$(awk '$4 == "CEN155" {count++} END {print count+0}' "$file")
    total=$((total + count))
    
    if [ $count -gt 0 ]; then
        echo "$sample: $count"
    fi
done

echo ""
echo "=== Execution Completed ==="
echo "Total CEN155 regions identified: $total"
echo "Output file: $OUTPUT_FILE"
echo ""
echo "Content preview (Top 10 rows):"
head -10 "$OUTPUT_FILE"
echo ""
echo "Statistics categorized by sample:"
tail -n +2 "$OUTPUT_FILE" | cut -f1 | sort | uniq -c | sort -rn