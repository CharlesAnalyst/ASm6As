library(limma)
library(edgeR)
# setwd("/home/galaxy/project/circRNA_analysis/find_circ/count_matrix")
# data <- read.table('count.txt',header=T,row.names = 1)
setwd("/home/galaxy/project/circRNA_analysis/sailfish-circ/count_matrix")
data <- read.table('TPM_matrix.txt',header=T,row.names = 2)
head(data)
#
dge <- DGEList(counts = data)
# apply scale normalization to read counts. TMM normalization method.
dge <- calcNormFactors(dge)
# convert the counts to logCPM to damp down the variances of logarithms of low counts.
logCPM <- cpm(dge, log = TRUE, prior.count = 3)
# Differential expression: limma-trend
group_list = c("wt", "ko", "wt", "ko")
design <- model.matrix(~group_list)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
trend_result <- topTable(fit, coef=ncol(design), number = Inf, sort.by = "P")
sum(trend_result$adj.P.Val<0.05)
write.table(trend_result, 'limma-trend_DE.txt',quote=F,sep="\t")
# Differential expression: voom
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
voom_result <- topTable(fit, coef = ncol(design), number = Inf, sort.by = "P")
sum(voom_result$adj.P.Val<0.05)
write.table(voom_result, 'limma-voom_DE.txt',quote=F,sep="\t")
# Differential expression: voom + sample-specific quality weights
# down-weight outlier samples
vwts <- voomWithQualityWeights(dge, design = design, normalization="none", plot = TRUE)
vfit2 <- lmFit(vwts)
vfit2 <- eBayes(vfit2)
combined_strategy_result <- topTable(vfit2, coef = 2, number = Inf, sort.by = "P")
sum(combined_strategy_result$adj.P.Val<0.05)
write.table(combined_strategy_result, 'limma-combined_strategy_DE.txt',quote=F,sep="\t")
#group_list = factor(c(rep('NT',3),rep('PA',3)))
# group_list = c("ko", "ko", "wt", "wt")
# design <- model.matrix(~group_list)
# colnames(design) <- levels(group_list)
# rownames(design) <- colnames(data)
# v <- voom(data,design,plot=TRUE,normalize="quantile")
# b <- v$E
# write.table(b,'limma_DE_result.txt',quote=F,sep="\t")
