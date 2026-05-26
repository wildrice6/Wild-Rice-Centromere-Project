#!/bin/bash
# Description: fixed_final_intersection.sh
# Revised Version: Utilizing CENH3 as the reference baseline to align CEN155 intervals.

CEN155_FILE="/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/all_CEN155_merged.tsv"
CENH3_FILE="/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/satellite_CENH3_region.tsv"
OUTPUT_DIR="/share/org/YZWL/yzwl_shahd/lxr/satellite_hor/fixed_final"

mkdir -p "$OUTPUT_DIR"

echo "=== Revised Analysis: Aligning CEN155 intervals based on CENH3 reference ==="
echo ""

# 1. Clean and standardize CENH3 data
echo "1. Cleaning and standardizing CENH3 data..."
CENH3_CLEAN="${OUTPUT_DIR}/cenh3_clean.tsv"

awk -F'\t' '
BEGIN {
    print "Sample\tChromosome\tStart\tEnd";
}
NR>1 {
    gsub(/["]/, "", $4);
    gsub(/["]/, "", $5);
    gsub(/,/, "", $4);
    gsub(/,/, "", $5);
    
    # Validate genomic interval coordinates
    if ($4 ~ /^[0-9]+$/ && $5 ~ /^[0-9]+$/ && $4 > 0 && $5 > $4) {
        printf "%s\t%s\t%s\t%s\n", $1, $3, $4, $5;
        valid_count++;
    } else {
        # Log invalid intervals for debugging
        invalid_count++;
        printf "Skipping invalid interval %d: %s %s %s-%s\n", invalid_count, $1, $3, $4, $5 > "/dev/stderr";
    }
}
END {
    printf "Valid intervals: %d, Invalid intervals skipped: %d\n", valid_count, invalid_count > "/dev/stderr";
}
' "$CENH3_FILE" > "$CENH3_CLEAN"

valid_cenh3_count=$(( $(wc -l < "$CENH3_CLEAN") - 1 ))
echo "  Total valid CENH3 intervals: $valid_cenh3_count"
echo ""

# 2. Process each CENH3 interval for intersection
echo "2. Processing individual CENH3 genomic intervals..."
FINAL_RESULT="${OUTPUT_DIR}/final_intersection.tsv"
DETAILED_RESULT="${OUTPUT_DIR}/detailed_intersection.tsv"
PROCESS_DIR="${OUTPUT_DIR}/process_files"
mkdir -p "$PROCESS_DIR"

echo -e "Sample\tChromosome\tStart\tEnd" > "$FINAL_RESULT"
echo -e "Sample\tChromosome\tCENH3_Start\tCENH3_End\tHas_CEN155\tFirst_CEN155_Start\tLast_CEN155_End\tIntersection_Start\tIntersection_End" > "$DETAILED_RESULT"

# Utilize temporary file to store cross-shell statistics
STATS_FILE="${OUTPUT_DIR}/temp_stats.txt"
echo "total_cenh3=0" > "$STATS_FILE"
echo "total_with_cen155=0" >> "$STATS_FILE"
echo "total_overlaps=0" >> "$STATS_FILE"

# Execute intersection logic via awk to maintain data integrity
awk -F'\t' -v cen155_file="$CEN155_FILE" -v process_dir="$PROCESS_DIR" \
    -v final_result="$FINAL_RESULT" -v detailed_result="$DETAILED_RESULT" '
# Load CENH3 standardized data (skip header)
FNR==NR && NR>1 {
    sample = $1;
    chr = $2;
    cenh3_start = $3;
    cenh3_end = $4;
    
    total_cenh3++;
    
    # Create intermediate processing log for each interval
    process_file = process_dir "/" sample "_" chr "_process.txt";
    printf "CENH3 Interval: %s %s %d-%d\n\n", sample, chr, cenh3_start, cenh3_end > process_file;
    
    # Query matching CEN155 intervals for the specific sample and chromosome
    cen155_cmd = "grep \"^" sample "\" \"" cen155_file "\" 2>/dev/null | grep \"\\t" chr "\\t\" 2>/dev/null | sort -k3,3n";
    cen155_cmd | getline cen155_line;
    close(cen155_cmd);
    
    # Retrieve count of available CEN155 intervals
    cen155_count_cmd = "grep \"^" sample "\" \"" cen155_file "\" 2>/dev/null | grep -c \"\\t" chr "\\t\" 2>/dev/null";
    cen155_count_cmd | getline cen155_count;
    close(cen155_count_cmd);
    
    printf "CEN155 regions identified: %d\n\n", cen155_count > process_file;
    
    if (cen155_count == 0) {
        # No CEN155 data available; utilize original CENH3 interval
        printf "No CEN155 data found; defaulting to CENH3 interval.\n" > process_file;
        printf "%s\t%s\t%s\t%s\n", sample, chr, cenh3_start, cenh3_end >> final_result;
        printf "%s\t%s\t%s\t%s\tNo\tNA\tNA\t%s\t%s\n", 
            sample, chr, cenh3_start, cenh3_end, cenh3_start, cenh3_end >> detailed_result;
        next;
    }
    
    total_with_cen155++;
    
    # Identify all CEN155 regions overlapping with the current CENH3 interval
    overlapping_count = 0;
    first_start = "";
    last_end = "";
    
    overlapping_file = process_dir "/" sample "_" chr "_overlapping.txt";
    printf "" > overlapping_file;
    
    # Iterate through and filter overlapping CEN155 regions
    while ((cen155_cmd | getline cen155_line) > 0) {
        split(cen155_line, cen155, "\t");
        cen155_start = cen155[3];
        cen155_end = cen155[4];
        
        # Intersection validation logic
        if (cen155_start <= cenh3_end && cen155_end >= cenh3_start) {
            overlapping_count++;
            total_overlaps++;
            
            # Calculate the boundaries of the intersection
            overlap_start = (cen155_start > cenh3_start) ? cen155_start : cenh3_start;
            overlap_end = (cen155_end < cenh3_end) ? cen155_end : cenh3_end;
            
            printf "%d. CEN155: %d-%d, Intersection part: %d-%d\n", 
                overlapping_count, cen155_start, cen155_end, overlap_start, overlap_end >> overlapping_file;
            
            # Update the global intersection boundaries for this interval
            if (first_start == "" || overlap_start < first_start) {
                first_start = overlap_start;
            }
            if (last_end == "" || overlap_end > last_end) {
                last_end = overlap_end;
            }
        }
    }
    close(cen155_cmd);
    
    printf "Number of overlapping CEN155 regions: %d\n\n", overlapping_count > process_file;
    
    if (overlapping_count > 0) {
        # Successful intersection: using the span from the first overlap start to the last overlap end
        printf "%d overlapping intervals confirmed.\n", overlapping_count > process_file;
        printf "Composite Intersection Start: %d\n", first_start > process_file;
        printf "Composite Intersection End: %d\n", last_end > process_file;
        printf "\n" > process_file;
        
        printf "Detailed overlapping interval data:\n" > process_file;
        system("cat \"" overlapping_file "\" >> \"" process_file "\"");
        printf "\n" > process_file;
        
        # Log results
        printf "%s\t%s\t%s\t%s\n", sample, chr, first_start, last_end >> final_result;
        printf "%s\t%s\t%s\t%s\tYes\t%s\t%s\t%s\t%s\n", 
            sample, chr, cenh3_start, cenh3_end, first_start, last_end, first_start, last_end >> detailed_result;
        
        # Log significant findings to console
        if (total_cenh3 <= 5 || overlapping_count > 0) {
            printf "  %s %s: Found %d overlapping CEN155 regions. Intersection: %d-%d\n", 
                sample, chr, overlapping_count, first_start, last_end > "/dev/stderr";
        }
    } else {
        # No overlaps found; utilize original CENH3 interval
        printf "No overlapping CEN155 intervals detected; utilizing original CENH3 interval.\n" > process_file;
        printf "%s\t%s\t%s\t%s\n", sample, chr, cenh3_start, cenh3_end >> final_result;
        printf "%s\t%s\t%s\t%s\tNo_Overlap\tNA\tNA\t%s\t%s\n", 
            sample, chr, cenh3_start, cenh3_end, cenh3_start, cenh3_end >> detailed_result;
        
        if (total_cenh3 <= 5) {
            printf "  %s %s: No overlaps found. Defaulting to CENH3: %d-%d\n", 
                sample, chr, cenh3_start, cenh3_end > "/dev/stderr";
        }
    }
    
    printf "Final Output Interval: %d-%d\n", 
        (first_start != "" ? first_start : cenh3_start), 
        (last_end != "" ? last_end : cenh3_end) > process_file;
    
    # Progress monitoring
    if (total_cenh3 % 10 == 0) {
        printf "Processed %d CENH3 intervals...\n", total_cenh3 > "/dev/stderr";
    }
}
END {
    # Store aggregated statistics
    printf "total_cenh3=%d\n", total_cenh3 > "'$STATS_FILE'";
    printf "total_with_cen155=%d\n", total_with_cen155 > "'$STATS_FILE'";
    printf "total_overlaps=%d\n", total_overlaps > "'$STATS_FILE'";
    
    printf "\nProcessing finalized:\n" > "/dev/stderr";
    printf "Total CENH3 intervals: %d\n", total_cenh3 > "/dev/stderr";
    printf "Intervals with corresponding CEN155 data: %d\n", total_with_cen155 > "/dev/stderr";
    printf "Total identified overlapping CEN155 regions: %d\n", total_overlaps > "/dev/stderr";
}
' "$CENH3_CLEAN"

# Retrieve statistics for the summary report
source "$STATS_FILE" 2>/dev/null

# 3. Generate summary report
echo ""
echo "3. Generating comprehensive analysis report..."
SUMMARY_FILE="${OUTPUT_DIR}/summary.txt"

# Calculate extended metrics
final_count=$(( $(wc -l < "$FINAL_RESULT") - 1 ))
detailed_count=$(( $(wc -l < "$DETAILED_RESULT") - 1 ))
has_overlap_count=$(tail -n +2 "$DETAILED_RESULT" | awk -F'\t' '$5=="Yes" {count++} END {print count+0}')
no_overlap_count=$(tail -n +2 "$DETAILED_RESULT" | awk -F'\t' '$5=="No_Overlap" {count++} END {print count+0}')
no_cen155_count=$(tail -n +2 "$DETAILED_RESULT" | awk -F'\t' '$5=="No" {count++} END {print count+0}')

cat > "$SUMMARY_FILE" << EOF
Analysis Report: CEN155 and CENH3 Genomic Intersection
======================================================
Report Generated: $(date)

Input Files:
- CEN155 Source: $(basename "$CEN155_FILE") ($(wc -l < "$CEN155_FILE") records)
- CENH3 Source: $(basename "$CENH3_FILE") ($(wc -l < "$CENH3_FILE") records)

Processing Statistics:
----------------------
Total CENH3 intervals analyzed: ${total_cenh3:-0}
Intervals with mapped CEN155 data: ${total_with_cen155:-0}
Intervals with confirmed CEN155 overlaps: $has_overlap_count
Total overlapping CEN155 sub-regions identified: ${total_overlaps:-0}

Output Files:
-------------
1. final_intersection.tsv - Consolidated intersection intervals (Primary Result)
   Format: Sample | Chromosome | Start | End
   Record Count: $final_count

2. detailed_intersection.tsv - Comprehensive intersection metrics
   Includes overlap status, CEN155 presence, and boundary details.
   Record Count: $detailed_count

3. summary.txt - This summary report.

4. process_files/ - Individual logs for each CENH3 interval processing step.
   Includes: 
   - {sample}_{chr}_process.txt - Step-by-step processing log
   - {sample}_{chr}_overlapping.txt - Specific overlapping regions (if applicable)

Categorization of Results:
--------------------------
1. Confirmed Overlaps: $has_overlap_count intervals
2. Presence of CEN155 but No Overlap: $no_overlap_count intervals
3. Absence of CEN155 Data: $no_cen155_count intervals

Sample Intervals with Identified Intersections (Top 10):
EOF

# Provide examples of confirmed intersections
if [ $has_overlap_count -gt 0 ]; then
    tail -n +2 "$DETAILED_RESULT" | awk -F'\t' '$5=="Yes" {print $1 " " $2 ": CENH3 " $3 "-" $4 ", Intersection " $8 "-" $9}' | head -10 >> "$SUMMARY_FILE"
else
    echo "No intersections identified." >> "$SUMMARY_FILE"
fi

echo "" >> "$SUMMARY_FILE"
echo "Primary Result Preview (Top 10):" >> "$SUMMARY_FILE"
head -10 "$FINAL_RESULT" >> "$SUMMARY_FILE"

# Cleanup temporary statistics file
rm -f "$STATS_FILE"

echo ""
echo "=== Execution Completed ==="
echo "Total CENH3 intervals analyzed: ${total_cenh3:-0}"
echo "Consolidated intersection intervals generated: $final_count"
echo "Detailed metrics records: $detailed_count"
echo ""
echo "Output Directory: $OUTPUT_DIR"
echo ""
echo "Key Output Files:"
ls -lh "$OUTPUT_DIR"/*.tsv "$OUTPUT_DIR"/*.txt 2>/dev/null | head -10
echo ""
echo "Total intermediate process logs generated: $(ls -la "$PROCESS_DIR"/* 2>/dev/null | wc -l)"
echo ""
echo "Final results preview (Top 10):"
head -10 "$FINAL_RESULT"