library(MeTPeak)
gtf <- "/data/database/GRCh38/GENCODE/gencode.v27.annotation.gtf"
ip_1 <- "/home/lzg/m6a_data/map_rate_60/SRR5573024.3.sorted.bam"
input_1 <- "/home/lzg/m6a_data/map_rate_60/SRR5573025.3.sorted.bam"
metpeak(GENE_ANNO_GTF=gtf,IP_BAM = ip_1,INPUT_BAM = input_1, EXPERIMENT_NAME="Test_sra")


library(exomePeak)
result <- exomepeak(GENE_ANNO_GTF = gtf, IP_BAM = c(ip_1), INPUT_BAM = c(input_1))

library(MeTPeak)
setwd('/home/galaxy/project/m6AQTL/data/RIP_seq/MeTPeak')
gtf <- '/data/database/GRCh38/GENCODE/gencode.v27.annotation.gtf'
filename<-read.table('/home/lzg/1905/Metpeak/input_ip_list.txt',sep = '\t',quote = NULL,stringsAsFactors=F)

#!/usr/bin/Script

for (i in nrow(filename) ){
  input<-filename[i,1]
  ip<-filename[i,2]
  IP_BAM <- unlist(strsplit(ip,','))
  INPUT_BAM <- unlist(strsplit(input,','))
  metpeak(GENE_ANNO_GTF=gtf,IP_BAM = IP_BAM,INPUT_BAM = INPUT_BAM,EXPERIMENT_NAME=paste(ip,input,sep = '_'))
}
