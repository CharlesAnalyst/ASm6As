library(edgeR)
setwd("/home/galaxy/project/circRNA_analysis/sailfish-circ/count_matrix")
data = read.table("reads_matrix.txt", sep = "\t", header = TRUE, row.names = "Name")
# df = data[, c("SRR1596086", "SRR1596088", "SRR1596098", "SRR1596100")]
# df = data[, c("heart_3_L6.bam", "heart_1_L6.bam", "heart_2_L7.bam", "brain_1_L6.bam", "brain_2_L7.bam", "brain_3_LX.bam")]
df = data
head(df)

group = c("WT", "KO", "WT", "KO")
cds <- DGEList(df, group = group)
y <- calcNormFactors(cds)
plotMDS(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)
# qlf.2vs1 <- glmQLFTest(fit, coef=2)
# topTags(qlf.2vs1)
qlf <- glmQLFTest(fit, coef = 1)
summary(decideTests(qlf))
topTags(qlf)
write.table(topTags(qlf, n = length(qlf)), sep = "\t", file="DEgenes_edgeR.txt", quote = FALSE)
