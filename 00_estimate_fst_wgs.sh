#!/bin/bash
## follow steps 01-06 in wgs/ to generate input vcf here
## author: sblain

#SBATCH --job-name=ref_concat
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --output=e01_ref_concat.%j

module load GCC/11.2.0 BCFtools/1.14

RUN_LOCATION="/scratch/user/sblain/get_AIMs/ref_panel/"

cd $RUN_LOCATION

#concatenate all files and only keep biallelic snps
bcftools concat ref_variantCalls/*ref.vcf.gz | \
	bcftools view -m2 -M2 -v snps | \
	bcftools filter -e 'F_MISSING > 0.25 || MAF <= 0.1' all_scaffolds.ref.vcf.gz \
		-o all_scaffolds.filtered.ref.vcf.gz

tabix all_scaffolds.filtered.ref.vcf.gz

echo "ref vcf concatenated and filtered"

module load VCFtools/0.1.16

echo "ref vcf filtered"

## estimate fst
vcftools --gzvcf all_scaffolds.filtered.ref.vcf.gz --weir-fst-pop inds_ref_inland --weir-fst-pop inds_ref_coastal --out ref

echo "ref fst calculated"

awk '{if($3>0.9 && $3!="-nan")print$0}' ref.weir.fst > ref.09.fst
awk '{if($3>0.2 && $3!="-nan")print$0}' ref.weir.fst > ref.02.fst
awk '{if($3>0.1 && $3!="-nan")print$0}' ref.weir.fst > ref.01.fst
