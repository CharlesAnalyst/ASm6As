library(Rtsne)
library(data.table)
library(preprocessCore)
library(sva)

setwd("/Data_Resource/m6A-seq/com_fetus_vs_adult/m6A_abundance/quantification/norm_by_scaleFactor")
m6a <- read.table("bedtools_quantification.matrix", row.names = 1)
m6a_t = t(m6a)
#test <- setDT(data.frame(m6a_t), keep.rownames = "Sample")[]
m6a_unique <- unique(m6a_t) # Remove duplicates
tsne_out <- Rtsne(as.matrix(m6a_unique), perplexity = 9) # Run TSNE
samples <- row.names(m6a_t)
# plot(tsne_out$Y, asp=1, type="p") # col=m6a_unique[,1],
x = tsne_out$Y[,1]
y = tsne_out$Y[,2]
group <- c("adult","adult","adult","adult","adult","adult","fetus", "fetus","fetus","adult","adult","fetus", "fetus","fetus","adult","fetus", "fetus","fetus","adult","adult","adult","adult","fetus","fetus","adult","adult","fetus","fetus")
pca_data <- data.frame(x, y, samples, group)
ggplot(pca_data, aes(x=x, y=y, colour=group))+geom_point()+geom_text(aes(label=samples),hjust=0, vjust=0)

########################## Normalize.quantiles ####################
dataMatNorm <- normalize.quantiles(m6a_t)
tsne_out <- Rtsne(as.matrix(dataMatNorm), perplexity = 9)
x = tsne_out$Y[,1]
y = tsne_out$Y[,2]
group <- c("adult","adult","adult","adult","adult","adult","fetus", "fetus","fetus","adult","adult","fetus", "fetus","fetus","adult","fetus", "fetus","fetus","adult","adult","adult","adult","fetus","fetus","adult","adult","fetus","fetus")
pca_data <- data.frame(x, y, samples, group)
ggplot(pca_data, aes(x=x, y=y, colour=group))+geom_point()+geom_text(aes(label=samples),hjust=0, vjust=0)


#######################################################################
######################### byGroup ####################################
#######################################################################
m6a <- read.table("bedtools_quantification_byGroup.matrix", row.names = 1)
m6a_t = t(m6a)
#test <- setDT(data.frame(m6a_t), keep.rownames = "Sample")[]
m6a_unique <- unique(m6a_t) # Remove duplicates
tsne_out <- Rtsne(as.matrix(m6a_unique), perplexity = 3) # Run TSNE
samples <- row.names(m6a_t)
# plot(tsne_out$Y, asp=1, type="p") # col=m6a_unique[,1],
x = tsne_out$Y[,1]
y = tsne_out$Y[,2]
group <- c("adult","fetus","adult", "fetus","adult","fetus","adult","fetus","adult","fetus")
pca_data <- data.frame(x, y, samples, group)
ggplot(pca_data, aes(x=x, y=y, colour=group))+geom_point()+geom_text(aes(label=samples),hjust=0, vjust=0)

############################# Normalize.quantiles ####################
dataMatNorm <- normalize.quantiles(m6a_t)
tsne_out <- Rtsne(as.matrix(dataMatNorm), perplexity = 3)
x = tsne_out$Y[,1]
y = tsne_out$Y[,2]
group <- c("adult","fetus","adult", "fetus","adult","fetus","adult","fetus","adult","fetus")
pca_data <- data.frame(x, y, samples, group)
ggplot(pca_data, aes(x=x, y=y, colour=group))+geom_point()+geom_text(aes(label=samples),hjust=0, vjust=0)


########################## Clustering #############################
setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASE_analysis/readCount/Fisher_test/sig/")
c <- read.table("brain_corr", row.names = 1, header = TRUE,na.strings = "NA", sep = "\t")
d <- as.dist(1-c)
hr <- hclust(d, method = "complete", members=NULL)
plot(hr, hang = -1)
plot(as.dendrogram(hr), edgePar=list(col=3, lwd=4), horiz=T)
