library(quantro)

setwd("/Data_Resource/m6A-seq/com_fetus_vs_adult/bams")

input_count <- read.table("input_samples.count", header = TRUE, row.names = 1)
head(input_count)
count <- as.matrix(input_count)
count <- log2(count)
types <- c("fetus","fetus", "fetus", "fetus", "fetus", "fetus", "fetus", "fetus", "fetus", "fetus", "fetus", "fetus", "fetus", "adult", "adult", "adult", "adult", "adult", "adult", "adult", "adult", "adult", "adult", "adult", "adult", "adult", "adult", "adult")
matdensity(count, groupFactor = types, xlab = " ", ylab = "density",main = "Raw reads count", brewer.n = 8, brewer.name = "Dark2")
#legend('right', c("adult","fetus"), col = c(1, 2), lty = 1, lwd = 3)
qtest <- quantro(object = count, groupFactor = types)

library(doParallel)
registerDoParallel(cores=20)
qtestPerm <- quantro(count, groupFactor = types, B = 1000)
qtestPerm

quantroPlot(qtestPerm)



####################################################################

setwd("/Data_Resource/m6A-seq/com_fetus_vs_adult/m6A_abundance/quantification/norm_by_scaleFactor")
input_count <- read.table("test.matrix", header = TRUE, row.names = 1)
head(input_count)

test <- input_count[, c(1:9)]

count <- as.matrix(input_count)
count <- log10(count)
types <- c("adult", "adult", "adult", "adult", "adult", "adult", "fetus","fetus", "fetus", "adult", "adult","fetus","fetus", "fetus","adult","fetus", "fetus","adult", "adult", "adult", "adult","fetus", "fetus", "adult", "adult","fetus", "fetus")
matdensity(count, groupFactor = types, xlab = " ", ylab = "density",main = "Raw reads count", brewer.n = 8, brewer.name = "Dark2")
#legend('right', c("adult","fetus"), col = c(1, 2), lty = 1, lwd = 3)
matboxplot(count, groupFactor = types, xlab = " ", ylab = "density",main = "Raw reads count", brewer.n = 8, brewer.name = "Dark2")

qtest <- quantro(object = count, groupFactor = types)

library(doParallel)
registerDoParallel(cores=20)
qtestPerm <- quantro(count, groupFactor = types, B = 1000)
qtestPerm

quantroPlot(qtestPerm)

################################################################

setwd("/home/galaxy/project/m6AQTL/2019_10_10/transcriptome_snp/SNPiR/shell/all_samples/results/quantification/")
input_count <- read.table("bedtools_quantification.matrix", header = TRUE, row.names = 1)
head(input_count)
count <- as.matrix(input_count)
count <- log10(count)
types <- rep("sample", 51)
matdensity(count, groupFactor = types, xlab = " ", ylab = "density",main = "RPKM ratio", brewer.n = 8, brewer.name = "Dark2")
#legend('right', c("adult","fetus"), col = c(1, 2), lty = 1, lwd = 3)
qtest <- quantro(object = count, groupFactor = types)

library(doParallel)
registerDoParallel(cores=20)
qtestPerm <- quantro(count, groupFactor = types, B = 1000)
qtestPerm

quantroPlot(qtestPerm)


