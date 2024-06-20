#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=144:0:0
#$ -l h_vmem=16G

module load samtools
module load bcftools

# File containing the list of URLs to download VCF files from
urls_file="/data/home/hmz251/WGS/REF/gnomadFilesChr.txt"

# Process each URL in the urls_file
while read -r url; do
    echo "Processing URL: $url"

    # Determine whether the file is for genomes or exomes
    if [[ "$url" == *"genomes"* ]]; then
        data_type="genomes"
    elif [[ "$url" == *"exomes"* ]]; then
        data_type="exomes"
    else
        echo "Data type not recognized."
        continue
    fi

    # Extract chromosome number from URL
    chr=$(echo "$url" | grep -oP 'chr[0-9XY]+')

    # Define filenames based on the extracted information
    filename="${data_type}.${chr}.vcf.bgz"
    output_vcf="filtered_${filename}"
    extracted_data="extracted_${filename%.vcf.bgz}.tsv"
    compressed_data="extracted_${filename%.vcf.bgz}.tsv.gz"

    # Download the VCF file using wget
    wget "$url" -O "$filename" -q

    # Filter VCF for AF != 0.00000 using bcftools
    bcftools view -i 'INFO/AF > 0.00000' "$filename" -o "$output_vcf"

    # Extract the desired fields: CHROM, POS, REF, ALT, and AF
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$output_vcf" > "$extracted_data"

    # Compress the extracted data
    gzip  "$extracted_data"

    # Remove the uncompressed extracted data, intermediate and original VCF files
    rm "$output_vcf" "$filename"

    echo "Finished processing $filename"
done < "$urls_file"

echo "All files processed."
