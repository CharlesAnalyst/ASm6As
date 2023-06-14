cutadapt -e 0.1 -n 1 -O 1 -q 10 -m 32 -a AGATCGGAAGAGCACACGTCT -A GATCGTCGGACTGTAGAACTC -o sample.R1.adapter3trim.fastq.gz -p sample.R2.adapter3trim.fastq.gz sample.R1.fastq.gz sample.R2.fastq.gz > sample.adapter3trim.log
zcat sample.R1.adapter3trim.fastq.gz | fastx_collapser -Q33 -i - -o sample.R1.adapter3trim.collapse.fa
zcat sample.R2.adapter3trim.fastq.gz | fastx_collapser -Q33 -i - -o sample.R2.adapter3trim.collapse.fa
cutadapt -e 0.1 -n 1 -O 1 -m 16 -u 5 -u -11 -o sample.R1.clean.fa sample.R1.adapter3trim.collapse.fa > sample.R1.barcoder5n3trim.log
cutadapt -e 0.1 -n 1 -O 1 -m 16 -u 11 -u -5 -o sample.R2.clean.fa sample.R2.adapter3trim.collapse.fa > sample.R2.barcoder5n3trim.log

STAR --genomeDir %s --twopassMode Basic --readFilesIn %s

gatk SplitNCigarReads -R hg38.fa -I bam -O split_bam 
gatk BaseRecalibrator -I split_bam -R hg38.fa -O recal_table --known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
gatk ApplyBQSR -I split_bam -R hg38.fa --output recal_bam -bqsr recal_table
gatk HaplotypeCaller -ERC GVCF -R hg38.fa -I recal_bam --dbsnp dbsnp_146.hg38.vcf.gz -O out.gvcf
gatk GenotypeGVCFs -R hg38.fa -V out.gvcf -O out.vcf
gatk VariantFiltration -R hg38.fa -V out.vcf -window 35 -cluster 3 --filter-name FS -filter "FS > 30.0" --filter-name QD -filter "QD < 2.0" -O filter.vcf
