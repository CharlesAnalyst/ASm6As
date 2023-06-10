library(pheatmap)
# setwd("/data5/galaxy/project/tf_analysis/motif_analysis/fimo_result")
# data <- read.table("zscore_filter-by-RPKM_result.txt", header = T, sep = "\t", row.names = 1)
setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/overlap_ASE/for_heatmap/")
# data <- read.table("Unsorted_diff_scale.txt", header = T, sep = "\t") # , row.names = 1
data <- read.table("Unsorted.txt2", header = T, sep = "\t")
col_list = c("ASm6A_ref", "ASm6A_alt") # "ASm6A_ref", "ASm6A_alt", "ASE_ref", "ASE_alt"
# col_list = c("ASE_ref", "ASE_alt")
# col_list = c("ASm6A_diff", "ASE_diff")
data <- data[col_list]
data = na.omit(data)
data <- log2(data+1)
head(data)
## gap row number from the python script: ASE_vs_ASm6A
pheatmap(data, show_rownames = FALSE, fontsize = 15, cluster_cols = FALSE, cluster_rows = FALSE,
         gaps_row = c(21029, 33563, 36490), filename = "ASm6A.pdf") # "ASm6A.pdf"
#########
# cutree_rows = 8,
# color = c("#F5F5F5", "#FA8072"), breaks = c(0, 1.5, 5)
cols <- makeColorRampPalette(c("white", "grey", "#ffdab9", "#F4A460"), 1.5/max(data), 100)
res <- pheatmap(data, cluster_cols = FALSE,  number_format = "%.4f", display_numbers = FALSE, show_rownames = FALSE, color = cols, cluster_rows = TRUE)
data.cluster <- cbind(data, cluster = cutree(res$tree_row, k = 7))
head(data.cluster)
write.table(data.cluster, file = "zscore_full_info.txt", sep = "\t", quote = FALSE)
###
file <- "/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/by_tissue/overlap_GWAS/enrichment/formatted_results.txt"
df <- read.table(file, sep = "\t", row.names = 1, header = TRUE)
head(df)
