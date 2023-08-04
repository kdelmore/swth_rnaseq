sbatch=$1
list=$2

while read prefix1	prefix2	prefix3
do
	echo "#!/bin/bash
#SBATCH --job-name=st.$prefix1.$prefix2
#SBATCH --ntasks=14
#SBATCH --ntasks-per-node=14
#SBATCH --time=6:00:00
#SBATCH --mem=360GB
#SBATCH --output=$prefix1.$prefix2.$prefix3.%j.out

module load GCC/9.3.0 VCFtools/0.1.16
bcftools='bcftools' ## not coded correctly in R code right now so can't change
mkdir tmp_"$prefix1"_"$prefix2"_"$prefix3"

## split into 5000000 bp regions

vcftools --gzvcf \$bcftools/bcftools_"$prefix1"_qual.vcf.gz \
--chr "$prefix1" \
--from-bp "$prefix2" \
--to-bp "$prefix3" \
--max-alleles 2 \
--recode --recode-INFO-all \
--out \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3"

## create pos file

awk -v OFS='\t' '{print\$1,\$2,\$4,\$5}' \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".recode.vcf | \
grep -v ^\#\# | \
sed '1d' > \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".pos.txt

## create gen file

vcftools --vcf \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".recode.vcf \
--keep ref_panel_drew.txt \
--012 \
--out \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3" \
--temp tmp_"$prefix1"_"$prefix2"_"$prefix3"

rm -rf tmp_"$prefix1"_"$prefix2"_"$prefix3"

awk -f trans.awk \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".012 > \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".gen.txt

sed -i 's/-1/NA/g' \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".gen.txt

sed -i '1 s/.*/KF14K01	KF15K01	KF15K02	LF11K02	LF09K01	JE28K04	JE28K05	JE28K06	JE28K07	JE30K05	JE30K06	JF07K01	JF07K02	JF08K03	JF08K05	JF08K08	KF14K02	KF26K02	KF26K05	KF26K07	KF26K09	KF26K13	KF26K14	KF27K01	KF27K04	KG06K04	KG06K06	KG06K07	KG07K02	KG07K03	KG07K04	KG09K01	LF16K01	LF06K01	AH25H03/' \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".gen.txt

sed -i 's/ /	/g' \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".gen.txt

module purge
module load Miniconda3/4.9.2
module load GCC/11.2.0 HTSlib/1.14
module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2

R -e 'library("STITCH");
options(bitmapType=\"cairo\");
STITCH(chr = \""$prefix1"\",
regionStart = "$prefix2",
regionEnd = "$prefix3",
buffer = 100000,
method = \"pseudoHaploid\",
switchModelIteration=37,
outputdir = \"imputed/"$prefix1"_"$prefix2"_"$prefix3"\",
posfile = \"bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".pos.txt\",
bamlist = \"bam_list.txt\",
genfile = \"bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3".gen.txt\",
K = 80,
tempdir = tempdir(),
nCores = 14,
nGen = 500,
shuffle_bin_radius=500,
iSizeUpperLimit=500000,
keepSampleReadsInRAM=TRUE,
outputSNPBlockSize=5000,
use_bx_tag=FALSE)'

module purge
module load GCC/11.2.0 BCFtools/1.14
tabix imputed/"$prefix1"_"$prefix2"_"$prefix3"/stitch."$prefix1"."$prefix2"."$prefix3".vcf.gz

rm \$bcftools/bcftools_"$prefix1"_"$prefix2"_"$prefix3"*
rm -rf imputed/"$prefix1"_"$prefix2"_"$prefix3"/input

" > $sbatch/$prefix1\_$prefix2\_$prefix3.sbatch

sbatch $sbatch/$prefix1\_$prefix2\_$prefix3.sbatch

done < $list
