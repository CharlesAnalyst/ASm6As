allelic_counts <- data.frame()
rtable_filenames <- list.files(path = "/home/galaxy/project/m6AQTL/2019_10_10/transcriptome_snp/SNPiR/shell/all_samples/results/ASE_analysis/ASE_readCount/input_count", pattern = "readcounts.txt", full.names = T)

for (filename in rtable_filenames){
  filename_base <- basename(filename);
  sample_name <- unlist(strsplit(filename_base, ".readcounts.txt"))[1];
  print(sample_name);
  allelic_counts_tmp <- read.table(filename, header = TRUE);
  allelic_counts_tmp$sample <-sample_name;
  allelic_counts <- rbind(allelic_counts, allelic_counts_tmp);
}

library(ggplot2)
ggplot(allelic_counts, aes(x=refCount, y=altCount)) + geom_point(alpha=0.1)

## remove homozygous SNPs (a lot of points with 0 ref or alt counts)
allelic_counts <- allelic_counts[allelic_counts$ref]
