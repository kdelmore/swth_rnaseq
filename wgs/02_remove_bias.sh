#!/bin/bash
# author: kira delmore
# date: feb 2020
# usage: ./remove_bias.sh <list>
# submits pbs file for each item in list

list="$1"
sbatch="$2"

while read prefix
do
        echo "!/bin/bash
#SBATCH --job-name="$prefix".bias
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=48000MB
#SBATCH --output=job.%j."$prefix".bias.out

## odds and ends to set

## directories

## directories
realign_bam_final_par1='/scratch/user/delmore/swth_low_coverage/bam_final_coastal'
realign_bam_final_par2='/scratch/user/delmore/swth_low_coverage/bam_final_inland'
path='/scratch/user/delmore/tools/'
log='/scratch/user/delmore/swth_low_coverage/bias_logs/'

## tools
module load GCC/8.2.0-2.31.1 SAMtools/1.9
module load Java/1.8.0_92
export JAVA_OPTS='-Xmx24g'

## remove reads that mapped at low quality and generate a list of those reads for each parent

samtools view -b -q30 \$realign_bam_final_par1/"$prefix".combo.bam > \$realign_bam_final_par1/"$prefix".combo.q30.bam 2> \$log/"$prefix".par1.q30.log
samtools view -F 4 \$realign_bam_final_par1/"$prefix".combo.q30.bam | cut -f 1 > \$realign_bam_final_par1/"$prefix"_pass.txt 2> \$log/"$prefix".par1.F4.log

samtools view -b -q30 \$realign_bam_final_par2/"$prefix".combo.bam > \$realign_bam_final_par2/"$prefix".combo.q30.bam 2> \$log/"$prefix".par2.q30.log
samtools view -F 4 \$realign_bam_final_par2/"$prefix".combo.q30.bam | cut -f 1 > \$realign_bam_final_par2/"$prefix"_pass.txt 2> \$log/"$prefix".par2.F4.log

## find intersection of reads

awk 'NR==FNR { lines[\$0]=1; next } \$0 in lines' \$realign_bam_final_par1/"$prefix"_pass.txt \$realign_bam_final_par2/"$prefix"_pass.txt > \$realign_bam_final_par1/"$prefix"_intersection.txt

## select those reads from the inland subspecies

\$path/ngsutilsj bam-filter --whitelist \$realign_bam_final_par1/"$prefix"_intersection.txt \$realign_bam_final_par2/"$prefix".combo.q30.bam \$realign_bam_final_par2/"$prefix".combo.q30.unbiased.bam
samtools index \$realign_bam_final_par2/"$prefix".combo.q30.unbiased.bam

## get final coverage statistics

samtools flagstat $realign_bam_final_par2/"$prefix".combo.q30.bam > $realign_bam_final_par2/"$prefix".combo.q30.flagstat.txt
samtools depth $realign_bam_final_par2/"$prefix".combo.q30.bam | awk '{sum+=\$3} END { print \"Average = \",sum/NR}' > $realign_bam_final_par2/"$prefix".combo.q30.depth.txt

samtools flagstat $realign_bam_final_par2/"$prefix".combo.unbiased.bam > $realign_bam_final_par2/"$prefix".combo.unbiased.flagstat.txt
samtools depth $realign_bam_final_par2/"$prefix".combo.unbiased.bam | awk '{sum+=\$3} END { print \"Average = \",sum/NR}' > $realign_bam_final_par2/"$prefix".combo.unbiased.depth.txt

" > $sbatch/$prefix.sh

sbatch $sbatch/$prefix.sh

done < $list
