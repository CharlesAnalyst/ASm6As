# 0 read cleaning and identification of allelic variants
After the read-cleaning and mapping steps according to the method([m6A-SAC-seq](https://github.com/shunliubio/m6A-SAC-seq)), the allelic variants were identified from the untreated data.

# 1 read cleaning and mapping of treated samples
After data cleaning, treated samples were mapped to the human genome according to [m6A-SAC-seq](https://github.com/shunliubio/m6A-SAC-seq).  Each sample could generate an alignment file.

# 2 assigned to reference and alternative allele
For each sample, two allelic bams were generated from the merged.bam file based on each allelic variants by asSeq.

# 3 get ASm6A sites
For each allelic variants, both allelic bam files were separately continued to run the pipeline([m6A-SAC-seq](https://github.com/shunliubio/m6A-SAC-seq)), and finally reports the single-base loci of m6A sites with corresponding mutation rate as modification level. Then a two tailed Fisher's exact test were carried out to assess the reads of m6A in REF and ALT files, and the significant allelic sites were considered as ASm6A sites.
