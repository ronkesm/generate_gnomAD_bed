#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Function to log messages with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Debug: Print all arguments
log_message "All arguments: $@"

# Check if enough arguments were passed
if [ "$#" -ne 3 ]; then
    log_message "Error: Incorrect number of arguments. Expected 3, got $#"
    log_message "Usage: $0 [VCF Directory] [AF Threshold or Filter File] [Temporary Filter Option]"
    exit 1
fi

# Directory containing VCF files
vcf_directory="$1"
# AF threshold or path to filter file
af_or_file="$2"
# Temporary filter option
use_temp_filter="$3"

# Debug: Print each argument
log_message "Argument 1 (VCF Directory): $vcf_directory"
log_message "Argument 2 (AF Threshold or Filter File): $af_or_file"
log_message "Argument 3 (Temporary Filter Option): $use_temp_filter"

# Path to the BED file (Adjust this as needed)
bed_file="/data/home/hmz251/WGS/BAM/VCF/temp_converted_bed_notation_corrected2.bed"

# Debug: Check if bed_file exists
if [ -f "$bed_file" ]; then
    log_message "BED file exists: $bed_file"
else
    log_message "Error: BED file does not exist: $bed_file"
    exit 1
fi

log_message "Starting script with parameters:"
log_message "VCF Directory: $vcf_directory"
log_message "AF Threshold or Filter File: $af_or_file"
log_message "Temporary Filter Option: $use_temp_filter"

# Debug: Check if vcf_directory exists
if [ -d "$vcf_directory" ]; then
    log_message "VCF directory exists: $vcf_directory"
else
    log_message "Error: VCF directory does not exist: $vcf_directory"
    exit 1
fi

# Determine whether af_or_file is a file or a numeric value
if [[ -f "$af_or_file" ]]; then
    filter_bed="$af_or_file"
    log_message "Using provided BED file for filtering: $filter_bed"
elif [[ "$af_or_file" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    filter_bed="temp_high_af_${af_or_file}.bed"
    log_message "Generating temporary filter for AF > $af_or_file"
    awk -v threshold="$af_or_file" '$6 > threshold' "$bed_file" > "$filter_bed"
    log_message "Generated temporary filter at $filter_bed"
    log_message "Number of entries in filter file: $(wc -l < "$filter_bed")"
else
    log_message "Error: Second argument must be a valid numeric AF threshold or an existing BED file path."
    exit 1
fi

# Debug: Check if filter_bed exists
if [ -f "$filter_bed" ]; then
    log_message "Filter BED file exists: $filter_bed"
else
    log_message "Error: Filter BED file does not exist: $filter_bed"
    exit 1
fi

# After the filter BED file check, add:
log_message "Starting to process VCF files"

# Process each VCF in the directory
vcf_count=0
log_message "Searching for VCF files in: $vcf_directory"
for vcf in "$vcf_directory"/*.vcf.gz; do
    log_message "Found VCF file: $vcf"
    if [[ -f "$vcf" ]]; then
        log_message "Processing file $vcf_count: $vcf"
        
        # Filter to only keep PASS variants
        pass_filtered_vcf="${vcf%.vcf.gz}.pass.vcf.gz"
        log_message "Filtering for PASS variants..."
        log_message "Running command: bcftools view -f PASS -O z -o \"$pass_filtered_vcf\" \"$vcf\""
        if bcftools view -f PASS -O z -o "$pass_filtered_vcf" "$vcf"; then
            log_message "bcftools command completed successfully"
        else
            log_message "Error: bcftools command failed"
            exit 1
        fi
        
        # Construct the new filename for the final output
        filtered_vcf="${vcf%.vcf.gz}.gnomad.filtered.vcf"
        
        log_message "Running bedtools intersect..."
        log_message "Command: bedtools intersect -a \"$pass_filtered_vcf\" -b \"$filter_bed\" -wa -wb"
        if bedtools intersect -b "$pass_filtered_vcf" -a "$filter_bed" -wa -wb > temp_intersect.bed; then
            log_message "bedtools intersect completed successfully"
        else
            log_message "Error: bedtools intersect failed"
            exit 1
        fi
        
        log_message "Running awk command..."
        if awk '$1==$7 && $2+1==$8 && $4==$10 && $5==$11 {print $7"\t"$8"\t"$8}' temp_intersect.bed > exact_matches.tmp; then
            log_message "awk command completed successfully"
        else
            log_message "Error: awk command failed"
            exit 1
        fi
        
        log_message "Number of exact matches: $(wc -l < exact_matches.tmp)"
        
        # Step 2: Use this list of exact matches to filter the original VCF
        log_message "Filtering VCF..."
        log_message "Command: bedtools intersect -v -a \"$pass_filtered_vcf\" -b exact_matches.tmp -header > \"$filtered_vcf\""
        if bedtools intersect -v -a "$pass_filtered_vcf" -b exact_matches.tmp -header > "$filtered_vcf"; then
            log_message "Final filtering completed successfully"
        else
            log_message "Error: Final filtering failed"
            exit 1
        fi
        
        log_message "Original VCF variant count: $(bcftools view -H "$pass_filtered_vcf" | wc -l)"
        log_message "Filtered VCF variant count: $(bcftools view -H "$filtered_vcf" | wc -l)"
        
        # Remove temporary files
        rm exact_matches.tmp temp_intersect.bed
        rm "$pass_filtered_vcf"
        
        log_message "Completed processing $vcf"
    else
        log_message "Warning: $vcf is not a file, skipping"
    fi
done

log_message "All VCF files processed. Total files processed: $vcf_count"

# End of script
log_message "Script completed successfully"