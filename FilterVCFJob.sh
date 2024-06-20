#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=144:0:0
#$ -l h_vmem=96G

module load bcftools
module load bedtools

bash FilterVCF_Gnomad_PASS.sh  /data/home/hmz251/WGS/BAM/VCF/strelkaBarts2 /data/home/hmz251/WGS/BAM/VCF/temp_high_af_0.0005.bed yes
