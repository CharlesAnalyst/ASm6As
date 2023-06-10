library(Trumpet)
setwd("/home/galaxy/project/m6AQTL/2019_10_10/transcriptome_snp/SNPiR/shell/leukemia_cell/results/samtools_sort/sample_rename")
f1 <- "leukemia_1.ip_leukemia.bam"
f2 <- "leukemia_2.ip_leukemia.bam"
f3 <- "leukemia_3.ip_leukemia.bam"
f4 <- "leukemia_4.ip_leukemia.bam"
f5 <- "leukemia_5.ip_leukemia.bam"

f6 <- "leukemia_1.input_leukemia.bam"
f7 <- "leukemia_2.input_leukemia.bam"
f8 <- "leukemia_3.input_leukemia.bam"
f9 <- "leukemia_4.input_leukemia.bam"
f10 <- "leukemia_5.input_leukemia.bam"

ip_bam <- c(f1,f3,f4,f5)
input_bam <- c(f6,f8,f9,f10)
contrast_ip_bam <- c(f2)
contrast_input_bam <- c(f7)
gtf <- "/data/database/GRCh38/GENCODE/gencode.v27.annotation.gtf"
trumpet_report <- Trumpet_report(IP_BAM = ip_bam, Input_BAM = input_bam, contrast_IP_BAM = contrast_ip_bam, contrast_Input_BAM = contrast_input_bam, condition1 = "untreated", condition2 = "treat", GENE_ANNO_GTF = gtf)
browseURL("Trumpet_report.html")


