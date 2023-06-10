library(biomaRt)
######################################
# get snp position according to rs id 
######################################
setwd("/home/galaxy/data/QTL/eQTL/sigle_tissue/GTEx_Analysis_v7_eQTL/rs_id")
result_dir = "/home/galaxy/data/QTL/eQTL/sigle_tissue/GTEx_Analysis_v7_eQTL/GRCh38_bed/"
## listEnsemblArchives()
# listMarts(archive=TRUE)
# listDatasets(mart = useMart(biomart = "ENSEMBL_MART_SNP"))
snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp") # Human Short Variants (SNPs and indels excluding flagged variants) (GRCh38.p12)
for (file in list.files()){
  print(file)
  get_snp_position(file)
}

get_snp_position <- function(in_file){
  df <- read.table(in_file, header = TRUE)
  head(df)
  snp_pos <- getBM(attributes = c("chr_name", "chrom_start", "chrom_end", "refsnp_id"),
                   filters = c("snp_filter"),
                   value = df$ID, 
                   mart = snpmart)
  result_file = paste(result_dir, in_file, sep = "")
  write.table(snp_pos, result_file, sep = "\t", quote = FALSE, row.names = FALSE)
}



# It's too slow!
######################################
# get rs id according to snp position
# use the previous version of GRCh37
######################################
setwd("/home/galaxy/data/QTL/eQTL/sigle_tissue/GTEx_Analysis_v7_eQTL/query_format")
result_dir = "/home/galaxy/data/QTL/eQTL/sigle_tissue/GTEx_Analysis_v7_eQTL/rs_id/"
# grch37 <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
# listMarts(grch37)
# listDatasets(grch37)
snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host="grch37.ensembl.org", path="/biomart/martservice")


for (file in list.files()){
  print(file)
  result_file = paste(result_dir, file, sep="")
  print(result_file)
  id_list <- c()
  df <- read.table(file, header = TRUE)
  # print(head(df))
  for(i in 1:nrow(df)){
    chr_name <- df[i, "chr"]
    start <- df[i, "start"]
    end <- df[i, "end"]
    raw_id <- df[i, "variant_id"]
    id <- get_snp_id(chr, start, end)
    i_result = paste(raw_id, ":", id)
    print(i_result)
  #   print(start)
  #   break
  # break
    id_list <- c(id_list, i_result)
  }
  lapply(id_list, write, result_file, append=TRUE)
}

get_snp_id <- function(chr, start, end){
  snp_id <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                  filters = c("chr_name", "start", "end"),
                  # values = list(1, 1161931, 1161932),
                  values = list(chr, start, end),
                  mart = snpmart)
  # print(snp_id)
  return(snp_id$refsnp_id)
}
