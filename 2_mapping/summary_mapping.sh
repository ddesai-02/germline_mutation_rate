#!/bin/bash

# Load project variables
source /home/devan/projects/rrg-shaferab/devan/germline_mutation_rate/variable.py

# Directory containing all *_summary.txt files
SUMMARY_DIR=$path/$sp/bam_files/summary

# Output file
OUTFILE=$SUMMARY_DIR/mapped_summary.txt

# Header
echo -e "sample\ttotal_reads\tproperly_paired\tpercent_mapped" > "$OUTFILE"

# Initialize totals for averaging
total_sum=0
mapped_sum=0
sample_count=0

# Process each summary file
for SUMFILE in "$SUMMARY_DIR"/*_summary.txt; do
    SAMPLE=$(basename "$SUMFILE" | sed 's/_summary.txt$//')

    TOTAL=$(grep "paired in sequencing" "$SUMFILE" | grep -Eo '^[0-9]+')
    MAPPED=$(grep "properly paired" "$SUMFILE" | grep -Eo '^[0-9]+')

    if [[ -z "$TOTAL" || -z "$MAPPED" ]]; then
        echo "Warning: Missing data in $SUMFILE. Skipping..."
        continue
    fi

    PERCENT=$(awk -v m=$MAPPED -v t=$TOTAL 'BEGIN { printf "%.2f", (m/t)*100 }')

    echo -e "${SAMPLE}\t${TOTAL}\t${MAPPED}\t${PERCENT}" >> "$OUTFILE"

    total_sum=$((total_sum + TOTAL))
    mapped_sum=$((mapped_sum + MAPPED))
    sample_count=$((sample_count + 1))
done

# Compute averages
if [ $sample_count -gt 0 ]; then
    avg_total=$(awk -v sum=$total_sum -v n=$sample_count 'BEGIN { printf "%.0f", sum/n }')
    avg_mapped=$(awk -v sum=$mapped_sum -v n=$sample_count 'BEGIN { printf "%.0f", sum/n }')
    avg_percent=$(awk -v m=$avg_mapped -v t=$avg_total 'BEGIN { printf "%.2f", (m/t)*100 }')

    echo -e "AVERAGE\t${avg_total}\t${avg_mapped}\t${avg_percent}" >> "$OUTFILE"
fi
