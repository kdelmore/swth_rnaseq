while read prefix
do

echo "SAMPLE = '"$prefix"'
BAM = 'ugh_bam_final_inland/"$prefix".combo.q30.unbiased.bam'
SEX = 'MALE'

FILE = open('bam_lists/regions_full_inland.txt')
for line in FILE:
	try:
		OUT = open('%s_%s_%s.lsf' %(SAMPLE, line.split(':')[0], line.split(':')[1].strip()), 'a')
	except IndexError:
		OUT = open('%s_%s.lsf' %(SAMPLE, line.strip()), 'a')
	OUT.write('#BSUB -L /bin/bash\n')
	OUT.write('#BSUB -J %s_%s\n' %(SAMPLE, line.strip()))
	OUT.write('#BSUB -n 1\n#BSUB -R "span[ptile=1]"\n#BSUB -R "rusage[mem=12700]"\n#BSUB -M 12700\n#BSUB -W 14:00\n')
	OUT.write('#BSUB -o stdout.%s_%s\n' %(SAMPLE, line.strip()))
	OUT.write('#BSUB -e stderr.%s_%s\n' %(SAMPLE, line.strip()))
	OUT.write('module load GATK/4.1.4.1-GCCcore-8.2.0-Python-3.7.2\n')
	if line.strip() == 'mt_inland':
		OUT.write('gatk HaplotypeCaller -I %s -O %s_%s.g.vcf.gz -R genomic_resources/inland_genome_folded.fasta -L %s --sample-name %s -ERC GVCF --sample-ploidy 1' %(BAM, BAM.strip('.bam'),line.strip(),  line.strip(), SAMPLE))
	elif line[0:4] == 'chrZ':
		if SEX == 'FEMALE':
			OUT.write('gatk HaplotypeCaller -I %s -O %s_%s.g.vcf.gz -R genomic_resources/inland_genome_folded.fasta -L %s --sample-name %s -ERC GVCF --sample-ploidy 1' %(BAM, BAM.strip('.bam'), line.strip(), line.strip(), SAMPLE))
		else:
			OUT.write('gatk HaplotypeCaller -I %s -O %s_%s.g.vcf.gz -R genomic_resources/inland_genome_folded.fasta -L %s --sample-name %s -ERC GVCF' %(BAM, BAM.strip('.bam'), line.strip(), line.strip(), SAMPLE))
	else:
		OUT.write('gatk HaplotypeCaller -I %s -O %s_%s.g.vcf.gz -R genomic_resources/inland_genome_folded.fasta -L %s --sample-name %s -ERC GVCF' %(BAM, BAM.strip('.bam'), line.strip(), line.strip(), SAMPLE))
	OUT.close()
" > hc_python_scripts/"$prefix".py

python hc_python_scripts/"$prefix".py
for i in "$prefix"_*.lsf; do bsub < $i; done
rm *"$prefix"_*.lsf
done < $1
