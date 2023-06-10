library(ggbiplot)
setwd("/home/galaxy/project/logistic_regression/overall_elevation/")
df <- read.table("expression_matrix_T.txt", header = TRUE, row.names = 1)
head(df)

pca <- prcomp(df, scale. = TRUE)
# str(pca)
ggbiplot(pca)

