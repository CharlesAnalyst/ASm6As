# 0 read cleaning and identification of allelic variants
After the read-cleaning and mapping steps according to the method([s41587-022-01487-9](https://www.nature.com/articles/s41587-022-01487-9)), the allelic variants were identified from the untreated data.

# 1 read cleaning and GLORI running
After data cleaning, GLORI treated samples were used to run GLORI analysis pipeline according to [GLORI-tools](https://github.com/liucongcas/GLORI-tools).  This step could generate an alignment file, merged.sorted.bam.

# 2 assigned to reference and alternative allele
Two allelic bams were generated from the merged.sorted.bam file based on each allelic variants by asSeq.

# 3 get ASm6A sites
For each allelic variants, both allelic bam files were separately continued to run the GLORI pipeline, and finally reports the single-base loci of m6A sites with corresponding A rate as modification level. Then a two tailed Fisher’s exact test were carried out to assess the reads of m6A in REF and ALT files, and the significant allelic sites were considered as ASm6A sites.
