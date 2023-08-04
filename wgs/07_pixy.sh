#!/bin/bash
#author: alex samano, 2022
#usage: dxy_pixy.sh <list>
#submits sbatch job for each scaffold in list
#see pixy documentation for formatting populations.txt

#SBATCH --job-name=dxy
#SBATCH --ntasks=2
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --mem=16GB
#SBATCH --error=job.dxy.err
#SBATCH --output=job.dxy.out

# filter sites

module load GCC/9.3.0 VCFtools/0.1.16

vcftools --gzvcf genotype_all_drew_allsites_2.vcf.gz \
--remove-indels \
--max-missing 0.8 \
--minDP 7 \
--maxDP 100 \
--recode --stdout | gzip -c > genotype_all_drew_allsites_2_filtered.vcf.gz

module load GCC/11.2.0 BCFtools/1.14

bgzip genotype_all_drew_allsites_2_filtered.vcf
tabix genotype_all_drew_allsites_2_filtered.vcf.gz

# use pixy to estimate parameters
module purge
module load Miniconda3/4.9.2
source activate pixy3

pixy --stats pi dxy fst \
--vcf genotype_all_drew_allsites_2_filtered.vcf.gz \
--populations populations.txt \
--bed_file genes.bed \
--output_prefix genes_dxy_fst

pixy --stats pi dxy fst \
--vcf genotype_all_drew_allsites_2_filtered.vcf.gz \
--populations populations.txt \
--bed_file genes_updown.bed \
--output_prefix genes_dxy_fst_updown
