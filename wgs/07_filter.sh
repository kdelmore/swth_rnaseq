#!/bin/bash
#SBATCH --job-name=filter_all
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32000MB
#SBATCH --output=job.%j.filter_all.out

## odds and ends
geno='/scratch/user/delmore/swth_ref_panel/haplotype_caller' ## directory
input='all' ## population
max_dep='100' ## +-3SD originally but strange coverage means i just put it to 100
max_miss='0.85' ## not sure

module load GCC/9.3.0 VCFtools/0.1.16
vcftools --gzvcf $geno/genotype_$input.vcf.gz \
--keep drew \
--recode --recode-INFO-all \
--out $geno/genotype_$input.drew 2> $geno/log_"$input"_drew

module purge
module load GCC/11.2.0 BCFtools/1.14
bcftools filter $geno/genotype_$input.drew.recode.vcf -g 5 --threads 6 --output $geno/genotype_$input.drew.indel5bp.vcf 2> $geno/log_"$input"_indel5bp

module purge
module load RTG-Tools/3.12.1-Java-1.8

rtg vcffilter -i $geno/genotype_$input.drew.indel5bp.vcf \
-o $geno/genotype_$input.drew.indel5bp.info.vcf \
--no-gzip \
--all-samples \
--keep-expr "INFO.QD > 10 && INFO.MQ > 40 && \
INFO.FS < 10 && INFO.SOR < 4 && INFO.ReadPosRankSum > -8 && \
INFO.MQRankSum > -12.5" --min-quality 20 2> $geno/log_"$input"_info

module purge
module load GCC/9.3.0 VCFtools/0.1.16

vcftools --vcf $geno/genotype_$input.drew.indel5bp.info.vcf \
--minDP 7 \
--maxDP "$max_dep" \
--recode --recode-INFO-all \
--out $geno/genotype_$input.drew.indel5bp.info.dp 2> $geno/log_"$input"_dp

module purge
module load RTG-Tools/3.12.1-Java-1.8

rtg vcffilter -i $geno/genotype_$input.drew.indel5bp.info.dp.recode.vcf \
-o $geno/genotype_$input.drew.indel5bp.info.dp.biallelic.vcf \
--snps-only \
--max-alleles 2 \
--all-samples --no-gzip 2> $geno/log_"$input"_biallelic

module purge
module load GCC/9.3.0 VCFtools/0.1.16

vcftools --vcf $geno/genotype_$input.drew.indel5bp.info.dp.biallelic.vcf \
--max-missing "$max_miss" \
--recode --recode-INFO-all \
--out $geno/genotype_$input.drew.indel5bp.info.dp.biallelic.missing 2> $geno/log_"$input"_missing
