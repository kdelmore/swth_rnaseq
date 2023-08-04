#!/bin/bash
# author: kira delmore
# date: jan 2019
# usage: ./01_trim_align.sh <list> <sbatch> <project> <ref>
# submits sbatch file for each item in list

## these things go on the command line 
list="$1" ## list of all the individuals
sbatch="$2" ## directory where you write the sbatch/lsf files
project="$3" ## this will be added to your bam files so use the same value for a single project
ref="$4" ## reference genome

while read prefix
do
        echo "#!/bin/bash
#SBATCH --job-name="$prefix".trim_align
#SBATCH --ntasks=2
#SBATCH --time=24:00:00
#SBATCH --mem=32000MB
#SBATCH --output=job.%j."$prefix".trim_align.out
 
## odds and ends to set
project='$3'
ref='$4'

## directories
TMPDIR=\$TMPDIR ## to read huge files that are temporary to a location other then your directory
raw='/scratch/user/delmore/swth_low_coverage/Delmore_04092021-245800561'
trim=\$TMPDIR
sam=\$TMPDIR
bam=\$TMPDIR
#log='/scratch/user/delmore/swth_low_coverage/trim_align_logs_inland/'
log='/scratch/user/delmore/swth_low_coverage/trim_align_logs_coastal/' ## you need to make these directories
#bam_final='/scratch/user/delmore/swth_low_coverage/ugh_bam_final_inland/'
bam_final='/scratch/user/delmore/swth_low_coverage/bam_final_coastal/'

## tools
module load GCCcore/9.3.0 Trim_Galore/0.6.6-Python-3.8.2

## trim reads
trim_galore --paired -fastqc --clip_R1 15 --clip_R2 15 --three_prime_clip_R1 5 --three_prime_clip_R2 5 --retain_unpaired \$raw/"$prefix"_R1.fastq \$raw/"$prefix"_R2.fastq -o \$trim

module purge
module load GCC/9.3.0  BWA/0.7.17

## align reads with bwa
bwa mem -M -t 8 \$ref \$trim/"$prefix2"_1_val_1.fq.gz \$trim/"$prefix2"_2_val_2.fq.gz > \$sam/"$prefix".sam 2> \$log/"$prefix".bwape.log 
bwa mem -M \$ref \$trim/"$prefix2"_1_unpaired_1.fq.gz > \$sam/"$prefix"_1_unpaired.sam 2> \$log/"$prefix".bwase1.log
bwa mem -M \$ref \$trim/"$prefix2"_2_unpaired_2.fq.gz > \$sam/"$prefix"_2_unpaired.sam 2> \$log/"$prefix".bwase2.log

module purge
module load GCC/8.2.0-2.31.1 SAMtools/1.9

## convert alignments from sam to bam
samtools view -Sb \$sam/"$prefix".sam > \$bam/"$prefix".bam 2> \$log/"$prefix".sampe.log
samtools view -Sb \$sam/"$prefix"_1_unpaired.sam > \$bam/"$prefix"_1_unpaired.bam 2> \$log/"$prefix".samse1.log
samtools view -Sb \$sam/"$prefix"_2_unpaired.sam > \$bam/"$prefix"_2_unpaired.bam 2> \$log/"$prefix".samse2.log

module purge
module load picard/2.18.27-Java-1.8

## clean up bam files
java -jar \$EBROOTPICARD/picard.jar CleanSam INPUT=\$bam/"$prefix".bam OUTPUT=\$bam/"$prefix".clean.bam  2> \$log/"$prefix".cleansam.log
java -jar \$EBROOTPICARD/picard.jar SortSam TMP_DIR=`pwd`/tmp INPUT=\$bam/"$prefix".clean.bam OUTPUT=\$bam/"$prefix".sort.bam MAX_RECORDS_IN_RAM=5000000 SORT_ORDER=coordinate 2> \$log/"$prefix".sortsam.log
java -jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups TMP_DIR=`pwd`/tmp I=\$bam/"$prefix".sort.bam O=\$bam/"$prefix".sortrg.bam MAX_RECORDS_IN_RAM=5000000 SORT_ORDER=coordinate RGID="$prefix" RGLB="$project" RGPL=ILLUMINA RGPU="$project" RGSM="$prefix" CREATE_INDEX=True 2> \$log/"$prefix".addRG.log
java -Xmx7g -jar \$EBROOTPICARD/picard.jar MarkDuplicates TMP_DIR=`pwd`/tmp INPUT=\$bam/"$prefix".sortrg.bam MAX_RECORDS_IN_RAM=5000000 OUTPUT=\$bam/"$prefix".duprem.bam M=\$log/"$prefix".duprem.log REMOVE_DUPLICATES=true 2> \$log/"$prefix".duprem.log

java -jar \$EBROOTPICARD/picard.jar CleanSam INPUT=\$bam/"$prefix"_1_unpaired.bam OUTPUT=\$bam/"$prefix"_1_unpaired.clean.bam 2> \$log/"$prefix"_1_unpaired.cleansam.log
java -jar \$EBROOTPICARD/picard.jar SortSam  INPUT=\$bam/"$prefix"_1_unpaired.clean.bam OUTPUT=\$bam/"$prefix"_1_unpaired.sort.bam SORT_ORDER=coordinate 2> \$log/"$prefix"_1_unpaired.sortsam.log
java -jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=\$bam/"$prefix"_1_unpaired.sort.bam O=\$bam/"$prefix"_1_unpaired.sortrg.bam SORT_ORDER=coordinate RGID="$prefix" RGLB="$project" RGPL=ILLUMINA RGPU="$project" RGSM="$prefix" CREATE_INDEX=True 2> \$log/"$prefix"_1_unpaired.addRG.log

java -jar \$EBROOTPICARD/picard.jar CleanSam INPUT=\$bam/"$prefix"_2_unpaired.bam OUTPUT=\$bam/"$prefix"_2_unpaired.clean.bam 2> \$log/"$prefix"_2_unpaired.cleansam.log
java -jar \$EBROOTPICARD/picard.jar SortSam  INPUT=\$bam/"$prefix"_2_unpaired.clean.bam OUTPUT=\$bam/"$prefix"_2_unpaired.sort.bam SORT_ORDER=coordinate 2> \$log/"$prefix"_2_unpaired.sortsam.log
java -jar \$EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=\$bam/"$prefix"_2_unpaired.sort.bam O=\$bam/"$prefix"_2_unpaired.sortrg.bam SORT_ORDER=coordinate RGID="$prefix" RGLB="$project" RGPL=ILLUMINA RGPU="$project" RGSM="$prefix" CREATE_INDEX=True 2> \$log/"$prefix"_2_unpaired.addRG.log

module purge
module load GCC/8.2.0-2.31.1 SAMtools/1.9

## generate final bams
samtools merge \$bam_final/"$prefix".combo.bam \$bam/"$prefix".duprem.bam \$bam/"$prefix"_1_unpaired.sortrg.bam \$bam/"$prefix"_2_unpaired.sortrg.bam > \$log/"$prefix".sammerge.log
samtools index \$bam_final/"$prefix".combo.bam

## check depth
samtools flagstat \$bam_final/"$prefix".combo.bam > \$bam_final/"$prefix".combo.flagstat.txt
samtools depth \$bam_final/"$prefix".combo.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > \$bam_final/"$prefix".combo.depth.txt

#echo \"Program finished with exit code \$? at: \`date\`\" >> \$JOBINFO

" > $sbatch/$prefix.sbatch

sbatch $sbatch/$prefix.sbatch

done < $list
