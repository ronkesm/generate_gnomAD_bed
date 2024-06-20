#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=144:0:0
#$ -l h_vmem=16G

# Define the directory containing the files
directory="/data/home/hmz251/WGS/REF"

# Define output BED file
output_bed="$directory/combined_bed_file.bed"

# Start by creating or clearing the existing output file
> "$output_bed"

# Process each gzipped TSV file for exomes
for file in "$directory"/extracted_exomes.chr*.tsv.gz; do
    echo "Processing $file for exomes..."
    zcat "$file" | awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2, $3, $4, $5}' >> "$output_bed"
done

# Process each gzipped TSV file for genomes
for file in "$directory"/extracted_genomes.chr*.tsv.gz; do
    echo "Processing $file for genomes..."
    zcat "$file" | awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2, $3, $4, $5}' >> "$output_bed"
done

echo "All files have been combined into $output_bed"