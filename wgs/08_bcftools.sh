#!/bin/bash
# author: kira delmore
# date: apr 2022
# usage: ./01_bcftools.sh <ref>

sbatch=$1
list=$2

while read prefix
do
        echo "#!/bin/bash
		
#SBATCH --job-name=bcftools
#SBATCH --output=job_bcftools_"$prefix"_%j
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load GCC/10.2.0 SAMtools/1.11 BCFtools/1.11

ref='/scratch/user/delmore/swth_ref/bCatUst1.pri.cur.20200806.MT.20190918_folded.fasta'
bcftools='bcftools'

bcftools mpileup \
-f \$ref  \
-b bam_list.txt \
--min-BQ 20 \
--min-MQ 20 \
--regions "$prefix" | \
bcftools call \
--skip-variants indels \
--multiallelic-caller \
--variants-only - > \$bcftools/bcftools_"$prefix".vcf

bgzip \$bcftools/bcftools_"$prefix".vcf
tabix \$bcftools/bcftools_"$prefix".vcf.gz

bcftools filter -O z -o \$bcftools/bcftools_"$prefix"_qual.vcf.gz -i '%QUAL>500' \$bcftools/bcftools_"$prefix".vcf.gz

" > $sbatch/$prefix.sbatch

sbatch $sbatch/$prefix.sbatch

done < $list
