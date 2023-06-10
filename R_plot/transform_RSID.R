library(biomaRt)
in_file = "/home/galaxy/project/m6AQTL/data/meqtl/total_data/abnormal_rs.txt"
snpmart = useMart(host="apr2019.archive.ensembl.org", biomart = "ENSEMBL_MART_SNP", dataset="hsapiens_snp")
data <- read.csv(in_file,sep = "\t", header = FALSE)
print(head(data))
for (row in 1:nrow(data)) {
  chr_name <- data[row, "V1"]
  pos  <- data[row, "V2"]
  print(chr_name)
  print(pos)
  print(getBM(attributes=c('refsnp_id','chrom_start'), filters = c('chr_name',"start","end"), values =list(chr_name,pos,pos),mart = snpmart))
}
