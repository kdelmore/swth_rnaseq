#!/bin/bash
#SBATCH --job-name=geno_all
#SBATCH --ntasks=1
#SBATCH --time=168:00:00
#SBATCH --mem=32000MB
#SBATCH --output=job.%j.geno_all.out

##tools

module load GCCcore/10.2.0 GATK/4.2.0.0-Java-11

## reference directories

ref='/scratch/user/delmore/swth_ref/bCatUst1.pri.cur.20200806.MT.20190918_folded.fasta'
db='/scratch/user/delmore/swth_ref_panel/all_db_drew/'
geno='/scratch/user/delmore/swth_ref_panel/haplotype_caller/'

## genotype GVCFS

gatk GenotypeGVCFs \
-R $ref \
-V gendb://$db \
-all-sites --max-alternate-alleles 4 --standard-min-confidence-threshold-for-calling 30 \
-O $geno/genotype_all_drew_allsites_2.vcf.gz
