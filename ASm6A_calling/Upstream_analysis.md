## SRR2fastq
```shell
SRRdir=

ls $SRRdir/*.1|while read id
do
        parallel-fastq-dump -t 10 -O ../fq --gzip -s ${id}
done
```

## trim
```shell

fq_dir=
trim_dir=

ls $fq_dir/*.gz|while read id
do
        /home/cjr/anaconda3/bin/trim_galore -q 20 --stringency 3 --length 30 $id --gzip -o $trim_dir
done
```

## alignment
```shell

ulimit -n 100000

prefix=$1

fq_dir=
fq=${fq_dir}/${prefix}.1_trimmed.fq.gz
genomeDir=
data_dir=

mkdir -p $data_dir
cd $data_dir

STAR --genomeDir $genomeDir --twopassMode Basic --readFilesIn $fq --runThreadN 20 --readFilesCommand zcat --outFileNamePrefix ./$prefix --outSAMtype BAM SortedByCoordinate --outSAMattrRGline ID:$prefix SM:$prefix PL:ILLUMINA LB:RNA PU:HiSeq2500
java -jar /data5/xlj/xlj/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${prefix}Aligned.sortedByCoord.out.bam O=${prefix}dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${prefix}output.metrics
samtools view -@ 10 -b -q 30 -o ${prefix}-uniq.bam ${prefix}dedupped.bam
samtools index ${prefix}-uniq.bam

```

## SNP calling
```shell
id=$1

bamdir=
hg19=hg19.fa
snpdb=1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
indeldb=Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
dbsnp=dbsnp_138.hg19.vcf.gz

snpdir=

outdir=${snpdir}/${id}
mkdir -p $outdir
bam=${bamdir}/${id}/${id}-uniq.bam
split_bam=${outdir}/${id}-split.bam
recal_table=${outdir}/${id}-recal.table
recal_bam=${outdir}/${id}-recal.bam
vcf=${outdir}/${id}_gatk.gvcf
filter_vcf=${outdir}/${id}_filter.vcf

/home/cjr/software/gatk-4.1.2.0/gatk SplitNCigarReads -R $hg19 -I $bam -O $split_bam
/home/cjr/software/gatk-4.1.2.0/gatk BaseRecalibrator -I $split_bam -R $hg19 -O $recal_table --known-sites $snpdb --known-sites $indeldb
/home/cjr/software/gatk-4.1.2.0/gatk ApplyBQSR -I $split_bam -R $hg19 --output $recal_bam -bqsr $recal_table
/home/cjr/software/gatk-4.1.2.0/gatk HaplotypeCaller -R $hg19 -I $recal_bam --dbsnp $dbsnp -O $vcf --dont-use-soft-clipped-bases true -stand-call-conf 20 --minimum-mapping-quality 20
/home/cjr/software/gatk-4.1.2.0/gatk VariantFiltration -R $hg19 -V $vcf -window 35 -cluster 3 --filter-name FS -filter "FS > 30.0" --filter-name QD -filter "QD < 2.0" -O $filter_vcf
```

## bam wasp
```shell
ulimit -n 10000

prefix=$1
sample=$2
ref_input=$3

genomeDir=

refSNP=${ref_input}/${ref_input}_filter.vcf
outdir=${prefix}/
inputfile=${prefix}.1_trimmed.fq.gz
mkdir -p $outdir
cd $outdir
STAR --genomeDir $genomeDir --twopassMode Basic --readFilesIn $inputfile --runThreadN 24 --readFilesCommand zcat --outFileNamePrefix ${outdir}/${prefix} --waspOutputMode SAMtag --outSAMtype BAM SortedByCoordinate --varVCFfile $refSNP --outSAMattrRGline ID:$prefix SM:$prefix

java -jar /data5/xlj/xlj/tools/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${prefix}Aligned.sortedByCoord.out.bam O=${prefix}dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${prefix}output.metric

samtools view -@ 10 -b -q 30 -o ${prefix}-uniq.bam ${prefix}dedupped.bam

samtools index ${prefix}-uniq.bam
```

## SNP filter
```python
#### remove repeat && RNA editing
snp_dir = ""
repeat = "RepeatMasker_hg19.record.bed"
editing = "Human_AG_all_hg19_v2.bed"
dbsnp = "dbsnp_146.hg19.bed"
####
file_l = os.listdir(snp_dir)
f_l = []
for f in file_l :
    if '.' in f :
        continue
    else :
        f_l.append(f)
f_l.sort()    
os.chdir(snp_dir)

for f in f_l :

    vcf = os.path.join(snp_dir, f, "%s_filter.vcf"%(f))
####
    result_dir = '%s/'%(f)
    result_file = os.path.join(result_dir, "%s_rm.vcf"%(f))
    os.system("bedtools intersect -a %s -b %s %s -v -wa -header > %s" % (vcf, repeat, editing, result_file))
    
#### overlapping with dbSNP
####
    result_file2 = os.path.join(result_dir, "%s_rm_dbsnp.vcf"%(f))
    os.system("bedtools intersect -a %s -b %s -wa -u -header > %s" % (result_file, dbsnp, result_file2))
```

## het and biallelic SNP
```shell
ls -d SRR* |while read id;do ls *rm_dbsnp.vcf | parallel 'cat {} | java -jar /home/software/snpEff/SnpSift.jar filter "(countHet() > 0)" > {.}_het.vcf';done

ls -d SRR* |while read id;do ls *rm_dbsnp_het.vcf | parallel "/home/cjr/software/gatk-4.1.2.0/gatk SelectVariants -R hg19.fa -V {} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O {.}_SNP.vcf"
```

## read counts
```shell
prefix = 
genome_fasta = 
$het_vcf_file = 

java -jar GenomeAnalysisTK.jar -R $genome_fasta -T ASEReadCounter -o ${prefix}.readcounts.txt -I bam_file -sites $het_vcf_file -U ALLOW_N_CIGAR_READS --minDepth 1 --minMappingQuality 255  --minBaseQuality 10 --allow_potentially_misencoded_quality_scores -L het_vcf_file
```
