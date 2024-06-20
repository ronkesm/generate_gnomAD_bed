# Generate gnomad bed file for filtering VCFs

## Why

Germline filtering is important for accurate SNV calling in whole genomes. Databases (eg. Gnomad) containing germline SNPs are routinely updated but not always provided in a user-friendly format. Gnomad files specifically are provided per-chromosome with a lot of extraneous information (creating very large files) unnecessary for routine SNV filtering. They are also split into SNPs identified from exomes versus genomes (non-overlapping). These scripts download the up to date (as of 20-Jun-2024 at least) gnomad v.4.1 database files for both exomes and genomes and generates a single bed file containing chromosome, ref + alt alleles, and the population allelic frequency. Example scripts to filter SNV files using the resultant bed file are also provided.

**Example use case:** The dndscv R package provides a sitednds function that computes dNdS ratios per specific mutation (rather than gene). We've found this function to be very sensitive to germline SNP contamination. In this case, filtering tumour VCFs using a more stringent population frequency cut-off allows us to produce a more refined analysis.

## How

This is a five-step process:
1. Download SNP VCFs from both exomes and genomes
2. Extract chr + alt + ref + AF and put into tsv files for each chromosome ( + remove any SNPs with an AF=0.00000) 
3. Remove the original files (they're huge. You may need more information from them - in which case - alter the `bcftools query` argument) 
4. Merge each tsv file containing the SNPs for each chromosome.
5. Filter VCFs with a custom population frequency cut-off (also filters out all non-`PASS` variants)

`GnomadV4.1.DownloadFilterAndExtractLoci.sh` takes care of steps 1-3, step 4 is done by `GnomadV4.1.MergeChromosomes.sh`, step 5 is done by `FilterVCF_Gnomad_PASS.sh` 



### Step-by-step guide
1. Provide gnomad URLs in `gnomadFilesChr.txt` (genomes) and `gnomadFilesChrExomes.txt` (exomes) eg:

```
# gnomadFilesChr.txt file
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr1.vcf.bgz
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr2.vcf.bgz
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr3.vcf.bgz
```

2. Run `GnomadV4.1.DownloadFilterAndExtractLoci.sh (and/or GnomadV4.1.DownloadFilterAndExtractLociExome)`. 
This script is formatted to be submitted to a HPC with 'qsub'. Replace the `urls_file` argument in the script with a path to the gnomadFilesChr(Exomes).txt file. This takes 24-48 hours.

3. Run `GnomadV4.1.MergeChromosomes.sh` in the same directory. 

4. You should now have a file called `CombinedExomeGenomeGnomadV4.1.Filter.bed` that looks like this (Size ~ 28.5 GB):

```
chr1    12137   12138   C       A       0.00909091
chr1    12194   12195   T       C       0.00374065
chr1    12197   12198   G       C       0.361991
chr1    12200   12201   C       G       0.00313808
chr1    12224   12225   C       T       0.000443262
```

6. Run `FilterVCF_Gnomad_PASS` providing a directory with VCFs, an AF threshold [OR] Bed file with SNP variants in the same format as `CombinedExomeGenomeGnomadV4.1.Filter.bed`. The last argument `(yes/no)` determines whether the new .bed file with your chosen AF threshold is removed after execution or not (so you can re-run the script without generating a new bed file). An example `qsub` script is provided in `FilterVCFJob.sh`. In your directory with VCFs, you should now have `[filename].filtered.vcf.gz` files

**Disclaimer: I used ChatGPT to help design _some_ of this code. You may want to verify that the script does what you want it to do before relying on its results.**
