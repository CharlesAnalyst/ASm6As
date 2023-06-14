trim_galore -q 20 --stringency 1 -e 0.3 --length 35
seqkit rmdup -s
fastx_trimmer -f 11
 
STAR --genomeDir %s --twopassMode Basic --readFilesIn %s

gatk SplitNCigarReads -R hg38.fa -I bam -O split_bam 
gatk BaseRecalibrator -I split_bam -R hg38.fa -O recal_table --known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
gatk ApplyBQSR -I split_bam -R hg38.fa --output recal_bam -bqsr recal_table
gatk HaplotypeCaller -ERC GVCF -R hg38.fa -I recal_bam --dbsnp dbsnp_146.hg38.vcf.gz -O out.gvcf
gatk GenotypeGVCFs -R hg38.fa -V out.gvcf -O out.vcf
gatk VariantFiltration -R hg38.fa -V out.vcf -window 35 -cluster 3 --filter-name FS -filter "FS > 30.0" --filter-name QD -filter "QD < 2.0" -O filter.vcf