library(effsize)
library(coin)
setwd("/home/galaxy/project/pausing_index/calculate_PI")
df = read.table("filter_data.txt", sep = "\t", header = TRUE)
GroupA <- df$log2WT
GroupB <- df$log2KO
result = cohen.d(GroupA, GroupB, pooled = TRUE, conf.level = 0.95)
print(result)
# cumulative distribution
plot(ecdf(GroupA), verticals = TRUE, pch = 46, col ="blue", xlab="Pausing Index", ylab="Probability", main="Empirical Cumulative Distribution")
lines(ecdf(GroupB), verticals = TRUE, pch = 46, col = "red")
legend(13, 0.8, legend = c("WT", "KO"), col = c("blue", "red"), lty = 1:1, cex = 0.8, box.lty=0)


library(tidyr)
library(ggplot2)
setwd("/home/galaxy/project/pausing_index/calculate_PI")
df = read.table("filter_data.txt", header = TRUE)
head(df)
df_long <- gather(df, condition, pausing_index, log2WT:log2KO, factor_key = TRUE)
head(df_long)
p0 <- ggplot(df_long, aes(x=condition, y =pausing_index)) + geom_boxplot()
ylim1 = boxplot.stats(df_long$pausing_index)$stats[c(1, 5)]
p1 = p0 + coord_cartesian(ylim = ylim1*1.05)
p1

