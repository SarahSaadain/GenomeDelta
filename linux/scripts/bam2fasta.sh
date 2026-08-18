#!/bin/bash
set -euo pipefail

# Check if correct number of arguments is provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <input_bam> <assembly> <min_cov> <min_len> <max_d> <output_folder> <prefix>"
    exit 1
fi

# Assign input arguments to variables
input_bam="$1"
assembly="$2"
min_cov="$3"
min_len="$4"
max_d="$5"
output_folder="$6"
prefix="$7"

# Create output folder if it doesn't exist
mkdir -p "$output_folder"

# Extract the file name from the prefix
filename="$prefix"

# Run samtools depth (gzipped, this file can get very large)
samtools depth -a "$input_bam" | gzip > "$output_folder/${filename}.bedgraph.gz"

# Filter for low coverage
zcat "$output_folder/${filename}.bedgraph.gz" | awk '$3 < '$min_cov > "$output_folder/${filename}-low_coverage.bedgraph"

# Create low coverage TSV
awk 'BEGIN{FS=OFS="\t"} {print $1, $2, $2}' "$output_folder/${filename}-low_coverage.bedgraph" > "$output_folder/${filename}-low_coverage.bed"

# Merge low coverage intervals
# First collapse touching/overlapping low-coverage bases into raw gaps
# (bedtools merge default distance is 0), then bridge raw gaps that are
# within max_d of each other. Counting how many raw gaps get bridged into
# each final region measures how "patchy" that gap is.
bedtools merge -i "$output_folder/${filename}-low_coverage.bed" > "$output_folder/${filename}-low_coverage_raw_gaps.bed"
bedtools merge -i "$output_folder/${filename}-low_coverage_raw_gaps.bed" -d "$max_d" -c 1 -o count > "$output_folder/${filename}-low_coverage_merged.tmp.bed"
# Remove one base from the start and end of each interval
awk 'BEGIN{OFS="\t"}{$2=$2+1; $3=$3-1}1' "$output_folder/${filename}-low_coverage_merged.tmp.bed" > "$output_folder/${filename}-low_coverage_merged.bed"
# Remove the temporary file
rm "$output_folder/${filename}-low_coverage_merged.tmp.bed"

# Remove small sequences (columns: chr, start, end, n_gaps_merged)
awk '($3 - $2) >= '$min_len "$output_folder/${filename}-low_coverage_merged.bed" > "$output_folder/${filename}-GD-with_gap_count.bed"
# credibility.sh/credibility.py expect a plain 3-column bed
cut -f1-3 "$output_folder/${filename}-GD-with_gap_count.bed" > "$output_folder/${filename}-GD.bed"

# Assign credibility to each sequence based on genomic context coverage
current_dir=$(dirname "$(readlink -f "$0")")
echo "Calculating coverage support for each gap"
bash "$current_dir/credibility.sh" "$output_folder/${filename}-GD.bed" "$output_folder/${filename}.bedgraph.gz" "$input_bam" "$assembly" "$output_folder" "$prefix"
python "$current_dir/credibility.py" "$output_folder/${filename}-GD.bed" "$output_folder/${filename}-GD-flanking.credibility" "$output_folder/${filename}.fai"

# Append the number of raw gaps merged into each region (patchiness) as the
# last column, and also bake it onto the name field (col 4) so it survives
# into the fasta header via `getfasta -name` and can be tracked per-cluster
# downstream (MSA2consensus.py, summarize_candidates.py).
paste "$output_folder/${filename}-GD-credibility.bed" <(cut -f4 "$output_folder/${filename}-GD-with_gap_count.bed") \
    | awk 'BEGIN{FS=OFS="\t"} {$4=$4"_"$6} 1' > "$output_folder/${filename}-GD-credibility.tmp.bed"
mv "$output_folder/${filename}-GD-credibility.tmp.bed" "$output_folder/${filename}-GD-credibility.bed"

# Use a loop to process each line in the BED file and get the corresponding cred_value
bedtools getfasta -fi "$assembly" -bed "$output_folder/${filename}-GD-credibility.bed" -fo "$output_folder/${filename}-GD.fasta" -name
sed -i'' -e 's/::.*$//' "$output_folder/${filename}-GD.fasta"

# Remove the temporary files
rm "$output_folder/${filename}-low_coverage.bedgraph"
rm "$output_folder/${filename}-low_coverage.bed"
rm "$output_folder/${filename}-low_coverage_raw_gaps.bed"
rm "$output_folder/${filename}-low_coverage_merged.bed"
rm "$output_folder/${filename}-GD-with_gap_count.bed"