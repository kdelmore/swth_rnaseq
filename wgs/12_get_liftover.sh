## chr names aren't the same between the chicken reference used in alignment and the bed file
sed 's/\r$//g' GCF_000002315.6_GRCg6a_assembly_report.txt | grep -v "^#" | cut -f1,5,7,10 > galGal_chr_key
 awk '{print$2,$1}' acckey > acckey2
./replace_chrs.pl acckey2 keep.bed > test.bed

## get coordinates in thrush
## halliftover is part of cactus
halLiftover --noDupes test.hal GRCg6a test.bed bCatUst1 bCatUst1.bed
