#!/bin/bash
#SBATCH --job-name=merge_db_all
#SBATCH --ntasks=1
#SBATCH --time=144:00:00
#SBATCH --mem=32000MB
#SBATCH --output=job.%j.db_all.out

module load GCCcore/10.2.0 GATK/4.2.0.0-Java-11
haplotype_caller='/scratch/user/delmore/swth_ref_panel/haplotype_caller'

gatk GenomicsDBImport \
-V $haplotype_caller/JE28K04.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JE28K05.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JE28K06.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JE28K07.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JE30K05.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JE30K06.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JF07K01.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JF07K02.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JF08K03.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JF08K05.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/JF08K08.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF14K01.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF14K02.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF15K01.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF15K02.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF26K02.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF26K05.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF26K07.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF26K09.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF26K13.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF26K14.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF27K01.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KF27K04.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KG06K04.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KG06K06.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KG06K07.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KG07K02.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KG07K03.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KG07K04.combo.q30.unbiased.g.vcf.gz \
-V $haplotype_caller/KG09K01.combo.q30.unbiased.g.vcf.gz \
-genomicsdb-workspace-path all_db_drew \
--intervals scaffold_100_arrow_ctg1 \
--intervals scaffold_101_arrow_ctg1 \
--intervals scaffold_102_arrow_ctg1 \
--intervals scaffold_103_arrow_ctg1 \
--intervals scaffold_104_arrow_ctg1 \
--intervals scaffold_105_arrow_ctg1 \
--intervals scaffold_106_arrow_ctg1 \
--intervals scaffold_107_arrow_ctg1 \
--intervals scaffold_108_arrow_ctg1 \
--intervals scaffold_109_arrow_ctg1 \
--intervals scaffold_10_arrow_ctg1:1-10000000 \
--intervals scaffold_10_arrow_ctg1:10000001-20000000 \
--intervals scaffold_10_arrow_ctg1:20000001-30000000 \
--intervals scaffold_10_arrow_ctg1:30000001-33537716 \
--intervals scaffold_111_arrow_ctg1 \
--intervals scaffold_112_arrow_ctg1 \
--intervals scaffold_113_arrow_ctg1 \
--intervals scaffold_114_arrow_ctg1 \
--intervals scaffold_115_arrow_ctg1 \
--intervals scaffold_116_arrow_ctg1 \
--intervals scaffold_117_arrow_ctg1 \
--intervals scaffold_118_arrow_ctg1 \
--intervals scaffold_119_arrow_ctg1 \
--intervals scaffold_11_arrow_ctg1:1-10000000 \
--intervals scaffold_11_arrow_ctg1:10000001-20000000 \
--intervals scaffold_11_arrow_ctg1:20000001-26454946 \
--intervals scaffold_120_arrow_ctg1 \
--intervals scaffold_121_arrow_ctg1 \
--intervals scaffold_122_arrow_ctg1 \
--intervals scaffold_123_arrow_ctg1 \
--intervals scaffold_124_arrow_ctg1 \
--intervals scaffold_125_arrow_ctg1 \
--intervals scaffold_126_arrow_ctg1 \
--intervals scaffold_127_arrow_ctg1 \
--intervals scaffold_128_arrow_ctg1 \
--intervals scaffold_129_arrow_ctg1 \
--intervals scaffold_12_arrow_ctg1:1-10000000 \
--intervals scaffold_12_arrow_ctg1:10000001-20000000 \
--intervals scaffold_12_arrow_ctg1:20000001-22310542 \
--intervals scaffold_130_arrow_ctg1 \
--intervals scaffold_131_arrow_ctg1 \
--intervals scaffold_132_arrow_ctg1 \
--intervals scaffold_133_arrow_ctg1 \
--intervals scaffold_134_arrow_ctg1 \
--intervals scaffold_135_arrow_ctg1 \
--intervals scaffold_136_arrow_ctg1 \
--intervals scaffold_137_arrow_ctg1 \
--intervals scaffold_138_arrow_ctg1 \
--intervals scaffold_139_arrow_ctg1 \
--intervals scaffold_13_arrow_ctg1:1-10000000 \
--intervals scaffold_13_arrow_ctg1:10000001-20000000 \
--intervals scaffold_13_arrow_ctg1:20000001-21745260 \
--intervals scaffold_140_arrow_ctg1 \
--intervals scaffold_141_arrow_ctg1 \
--intervals scaffold_142_arrow_ctg1 \
--intervals scaffold_143_arrow_ctg1 \
--intervals scaffold_144_arrow_ctg1 \
--intervals scaffold_145_arrow_ctg1 \
--intervals scaffold_146_arrow_ctg1 \
--intervals scaffold_147_arrow_ctg1 \
--intervals scaffold_148_arrow_ctg1 \
--intervals scaffold_149_arrow_ctg1 \
--intervals scaffold_14_arrow_ctg1:1-10000000 \
--intervals scaffold_14_arrow_ctg1:10000001-20000000 \
--intervals scaffold_14_arrow_ctg1:20000001-20997689 \
--intervals scaffold_150_arrow_ctg1 \
--intervals scaffold_151_arrow_ctg1 \
--intervals scaffold_152_arrow_ctg1 \
--intervals scaffold_153_arrow_ctg1 \
--intervals scaffold_154_arrow_ctg1 \
--intervals scaffold_155_arrow_ctg1 \
--intervals scaffold_156_arrow_ctg1 \
--intervals scaffold_157_arrow_ctg1 \
--intervals scaffold_158_arrow_ctg1 \
--intervals scaffold_159_arrow_ctg1 \
--intervals scaffold_160_arrow_ctg1 \
--intervals scaffold_161_arrow_ctg1 \
--intervals scaffold_162_arrow_ctg1 \
--intervals scaffold_163_arrow_ctg1 \
--intervals scaffold_164_arrow_ctg1 \
--intervals scaffold_165_arrow_ctg1 \
--intervals scaffold_166_arrow_ctg1 \
--intervals scaffold_167_arrow_ctg1 \
--intervals scaffold_168_arrow_ctg1 \
--intervals scaffold_169_arrow_ctg1 \
--intervals scaffold_16_arrow_ctg1:1-10000000 \
--intervals scaffold_16_arrow_ctg1:10000001-19392598 \
--intervals scaffold_170_arrow_ctg1 \
--intervals scaffold_171_arrow_ctg1 \
--intervals scaffold_172_arrow_ctg1 \
--intervals scaffold_173_arrow_ctg1 \
--intervals scaffold_174_arrow_ctg1 \
--intervals scaffold_175_arrow_ctg1 \
--intervals scaffold_176_arrow_ctg1 \
--intervals scaffold_177_arrow_ctg1 \
--intervals scaffold_178_arrow_ctg1 \
--intervals scaffold_179_arrow_ctg1 \
--intervals super_scaffold_17:1-10000000 \
--intervals super_scaffold_17:10000001-20000000 \
--intervals super_scaffold_17:20000001-21568194 \
--intervals scaffold_180_arrow_ctg1 \
--intervals scaffold_181_arrow_ctg1 \
--intervals scaffold_182_arrow_ctg1 \
--intervals scaffold_183_arrow_ctg1 \
--intervals scaffold_184_arrow_ctg1 \
--intervals scaffold_185_arrow_ctg1 \
--intervals scaffold_186_arrow_ctg1 \
--intervals scaffold_187_arrow_ctg1 \
--intervals scaffold_188_arrow_ctg1 \
--intervals scaffold_189_arrow_ctg1 \
--intervals scaffold_18_arrow_ctg1:1-10000000 \
--intervals scaffold_18_arrow_ctg1:10000001-16676293 \
--intervals scaffold_190_arrow_ctg1 \
--intervals scaffold_191_arrow_ctg1 \
--intervals scaffold_192_arrow_ctg1 \
--intervals scaffold_19_arrow_ctg1:1-10000000 \
--intervals scaffold_19_arrow_ctg1:10000001-15618554 \
--intervals super_scaffold_1:1-10000000 \
--intervals super_scaffold_1:10000001-20000000 \
--intervals super_scaffold_1:20000001-30000000 \
--intervals super_scaffold_1:30000001-40000000 \
--intervals super_scaffold_1:40000001-50000000 \
--intervals super_scaffold_1:50000001-60000000 \
--intervals super_scaffold_1:60000001-70000000 \
--intervals super_scaffold_1:70000001-80000000 \
--intervals super_scaffold_1:80000001-90000000 \
--intervals super_scaffold_1:90000001-100000000 \
--intervals super_scaffold_1:100000001-110000000 \
--intervals super_scaffold_1:110000001-120000000 \
--intervals super_scaffold_1:120000001-130000000 \
--intervals super_scaffold_1:130000001-140000000 \
--intervals super_scaffold_1:140000001-150000000 \
--intervals super_scaffold_1:150000001-165774257 \
--intervals scaffold_20_arrow_ctg1:1-10000000 \
--intervals scaffold_20_arrow_ctg1:10000001-14366659 \
--intervals scaffold_22_arrow_ctg1:1-10000000 \
--intervals scaffold_22_arrow_ctg1:10000001-12176784 \
--intervals scaffold_23_arrow_ctg1:1-10000000 \
--intervals scaffold_23_arrow_ctg1:10000001-11336639 \
--intervals scaffold_24_arrow_ctg1:1-10000000 \
--intervals scaffold_24_arrow_ctg1:10000001-11098648 \
--intervals super_scaffold_25:1-10000000 \
--intervals super_scaffold_25:10000001-12285394 \
--intervals scaffold_26_arrow_ctg1 \
--intervals scaffold_27_arrow_ctg1 \
--intervals scaffold_28_arrow_ctg1 \
--intervals scaffold_29_arrow_ctg1 \
--intervals scaffold_2_arrow_ctg1:1-10000000 \
--intervals scaffold_2_arrow_ctg1:10000001-20000000 \
--intervals scaffold_2_arrow_ctg1:20000001-30000000 \
--intervals scaffold_2_arrow_ctg1:30000001-40000000 \
--intervals scaffold_2_arrow_ctg1:40000001-50000000 \
--intervals scaffold_2_arrow_ctg1:50000001-60000000 \
--intervals scaffold_2_arrow_ctg1:60000001-70000000 \
--intervals scaffold_2_arrow_ctg1:70000001-80000000 \
--intervals scaffold_2_arrow_ctg1:80000001-90000000 \
--intervals scaffold_2_arrow_ctg1:90000001-100000000 \
--intervals scaffold_2_arrow_ctg1:100000001-110000000 \
--intervals scaffold_2_arrow_ctg1:110000001-120000000 \
--intervals scaffold_2_arrow_ctg1:120000001-124776741 \
--intervals scaffold_30_arrow_ctg1 \
--intervals scaffold_31_arrow_ctg1 \
--intervals scaffold_32_arrow_ctg1 \
--intervals super_scaffold_33 \
--intervals super_scaffold_35 \
--intervals super_scaffold_38 \
--intervals super_scaffold_39 \
--intervals super_scaffold_3:1-10000000 \
--intervals super_scaffold_3:10000001-20000000 \
--intervals super_scaffold_3:20000001-30000000 \
--intervals super_scaffold_3:30000001-40000000 \
--intervals super_scaffold_3:40000001-50000000 \
--intervals super_scaffold_3:50000001-60000000 \
--intervals super_scaffold_3:60000001-70000000 \
--intervals super_scaffold_3:70000001-80000000 \
--intervals super_scaffold_3:80000001-90000000 \
--intervals super_scaffold_3:90000001-100000000 \
--intervals super_scaffold_3:100000001-110000000 \
--intervals super_scaffold_3:110000001-120000000 \
--intervals super_scaffold_3:120000001-120187523 \
--intervals super_scaffold_45 \
--intervals scaffold_4_arrow_ctg1:1-10000000 \
--intervals scaffold_4_arrow_ctg1:10000001-20000000 \
--intervals scaffold_4_arrow_ctg1:20000001-30000000 \
--intervals scaffold_4_arrow_ctg1:30000001-40000000 \
--intervals scaffold_4_arrow_ctg1:40000001-50000000 \
--intervals scaffold_4_arrow_ctg1:50000001-60000000 \
--intervals scaffold_4_arrow_ctg1:60000001-70000000 \
--intervals scaffold_4_arrow_ctg1:70000001-76913702 \
--intervals scaffold_50_arrow_ctg1 \
--intervals scaffold_51_arrow_ctg1 \
--intervals scaffold_54_arrow_ctg1 \
--intervals super_scaffold_53 \
--intervals super_scaffold_55 \
--intervals super_scaffold_56 \
--intervals super_scaffold_58 \
--intervals super_scaffold_59 \
--intervals super_scaffold_5:1-10000000 \
--intervals super_scaffold_5:10000001-20000000 \
--intervals super_scaffold_5:20000001-30000000 \
--intervals super_scaffold_5:30000001-40000000 \
--intervals super_scaffold_5:40000001-50000000 \
--intervals super_scaffold_5:50000001-60000000 \
--intervals super_scaffold_5:60000001-70000000 \
--intervals super_scaffold_5:70000001-77026709 \
--intervals scaffold_62_arrow_ctg1 \
--intervals scaffold_68_arrow_ctg1 \
--intervals super_scaffold_46 \
--intervals scaffold_6_arrow_ctg1:1-10000000 \
--intervals scaffold_6_arrow_ctg1:10000001-20000000 \
--intervals scaffold_6_arrow_ctg1:20000001-30000000 \
--intervals scaffold_6_arrow_ctg1:30000001-40000000 \
--intervals scaffold_6_arrow_ctg1:40000001-50000000 \
--intervals scaffold_6_arrow_ctg1:50000001-60000000 \
--intervals scaffold_6_arrow_ctg1:60000001-63528738 \
--intervals scaffold_72_arrow_ctg1 \
--intervals scaffold_73_arrow_ctg1 \
--intervals scaffold_74_arrow_ctg1 \
--intervals scaffold_75_arrow_ctg1 \
--intervals scaffold_78_arrow_ctg1 \
--intervals scaffold_79_arrow_ctg1 \
--intervals super_scaffold_7:1-10000000 \
--intervals super_scaffold_7:10000001-20000000 \
--intervals super_scaffold_7:20000001-30000000 \
--intervals super_scaffold_7:30000001-40000000 \
--intervals super_scaffold_7:40000001-50000000 \
--intervals super_scaffold_7:50000001-60000000 \
--intervals super_scaffold_7:60000001-68462074 \
--intervals scaffold_80_arrow_ctg1 \
--intervals scaffold_84_arrow_ctg1 \
--intervals scaffold_85_arrow_ctg1 \
--intervals scaffold_86_arrow_ctg1 \
--intervals scaffold_87_arrow_ctg1 \
--intervals scaffold_89_arrow_ctg1 \
--intervals scaffold_8_arrow_ctg1:1-10000000 \
--intervals scaffold_8_arrow_ctg1:10000001-20000000 \
--intervals scaffold_8_arrow_ctg1:20000001-30000000 \
--intervals scaffold_8_arrow_ctg1:30000001-40000000 \
--intervals scaffold_8_arrow_ctg1:40000001-40436149 \
--intervals scaffold_90_arrow_ctg1 \
--intervals scaffold_91_arrow_ctg1 \
--intervals scaffold_92_arrow_ctg1 \
--intervals scaffold_93_arrow_ctg1 \
--intervals scaffold_94_arrow_ctg1 \
--intervals scaffold_95_arrow_ctg1 \
--intervals scaffold_96_arrow_ctg1 \
--intervals scaffold_97_arrow_ctg1 \
--intervals scaffold_98_arrow_ctg1 \
--intervals scaffold_99_arrow_ctg1 \
--intervals scaffold_9_arrow_ctg1:1-10000000 \
--intervals scaffold_9_arrow_ctg1:10000001-20000000 \
--intervals scaffold_9_arrow_ctg1:20000001-30000000 \
--intervals scaffold_9_arrow_ctg1:30000001-36548054 \
--intervals mtDNA
