#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=144:0:0
#$ -l h_vmem=16G

# Check if enough arguments were passed
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 [VCF Directory] [AF Threshold or Filter File] [Temporary Filter Option]"
    exit 1
fi

# Directory containing VCF files
vcf_directory="$1"

# AF threshold or path to filter file
af_or_file="$2"

# Temporary filter option
use_temp_filter="$3"

# Path to the BED file (Adjust this as needed)
bed_file="/data/home/hmz251/WGS/REF/CombinedExomeGenomeGnomadV4.1.Filter.bed"


# Determine whether af_or_file is a file or a numeric value
if [[ -f "$af_or_file" ]]; then
    filter_bed="$af_or_file"
    echo "Using provided BED file for filtering: $filter_bed"
elif [[ "$af_or_file" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    # Create a temporary BED file with entries having AF > provided threshold
    # Naming the temporary file with the AF threshold used for easy identification
    filter_bed="temp_high_af_${af_or_file}.bed"
    awk -v threshold="$af_or_file" '$6 > threshold' "$bed_file" > "$filter_bed"
    echo "Generated temporary filter for AF > $af_or_file at $filter_bed"
else
    echo "Error: Second argument must be a valid numeric AF threshold or an existing BED file path."
    exit 1
fi

# Process each VCF in the directory
for vcf in "$vcf_directory"/*.vcf.gz; do
    if [[ -f "$vcf" ]]; then  # Check if it's a file
        echo "Processing $vcf..."

        # Filter to only keep PASS variants
        pass_filtered_vcf="${vcf%.vcf.gz}.pass.vcf.gz"
        bcftools view -f PASS -O z -o "$pass_filtered_vcf" "$vcf"

        # Construct the new filename for the final output
        filtered_vcf="${vcf%.vcf.gz}.gnomad.filtered.vcf"

        # Use bedtools to exclude variants present in the filter BED file
        echo "Running bedtools command:"
        echo "bedtools intersect -v -a \"$pass_filtered_vcf\" -b \"$filter_bed\" -header > \"$filtered_vcf\""
   


        bedtools intersect -v -a "$pass_filtered_vcf" -b "$filter_bed" -header > "$filtered_vcf" 2>intersect_error.log
        if [[ $? -ne 0 ]]; then
            echo "bedtools intersect failed. Check intersect_error.log for details."
            cat intersect_error.log
            exit 1
        else
            echo "bedtools intersect succeeded."
            ls -lh "$filtered_vcf"  # Check the output file size and permissions
        fi

        # Optionally, remove the intermediate PASS-filtered file
        rm "$pass_filtered_vcf"
    fi
done

# Optionally, remove the temporary BED file if not specified to keep
if [[ "$use_temp_filter" == "no" && -f "$filter_bed" && "$filter_bed" == temp_high_af* ]]; then
    rm "$filter_bed"
    echo "Temporary filter removed: $filter_bed"
fi

echo "All VCF files processed."
