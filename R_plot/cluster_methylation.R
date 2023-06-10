library(tibble)
library(corrplot)
base_dir = "/home/galaxy/project/m6AQTL/2019_10_10/transcriptome_snp/SNPiR/shell/HeLa_cell/combined_results/results/MatrixEQTL_result/"
in_file = paste(base_dir, "hg19_quantification_T_R.matrix", sep = "")
df = read.csv(in_file)

M = cor(df)

mydata <- na.omit(df)
mydata <- scale(mydata)
fit <- kmeans(mydata, 2)
d <- dist(mydata, method = "euclidean")
fit <- hclust(d, method = "ward")
plot(fit)

#library(ggfortify)
# https://cmdlinetips.com/2019/04/introduction-to-pca-with-r-using-prcomp/
###
setwd("/Data_Resource/m6A-seq/com_fetus_vs_adult/m6A_abundance/quantification")
df = read.table("bedtools_quantification.matrix2", sep = "\t", row.names = 1, header = T)
pca= prcomp(log2(df+1) , center=T, scale=T)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
xlab = paste0("PC1: ", round(var_explained[1]*100,1), "% variance")
ylab = paste0("PC2: ", round(var_explained[2]*100,1), "% variance")
plot(pca$rotation[,1],pca$rotation[,2], xlab = xlab, ylab = ylab, labels=rownames(pca)) + text(pca$rotation[,1], pca$rotation[,2], labels=lables, cex= 0.7, pos = 3)

#
# pca_res$x %>% 
# as.data.frame %>%
#   ggplot(aes(x=PC1,y=PC2)) + geom_point(size=4) +
#   theme_bw(base_size=32) + 
#   labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#        y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#   theme(legend.position="top")
#