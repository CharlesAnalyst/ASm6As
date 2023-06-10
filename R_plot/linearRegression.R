###  RBPs for which changes in binding affinities were significantlly correlated with ASm6A imbalance
setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/motifbreakR/by_ind/plot/")
#### SRSF9
x=c(0.817805, 0.943519, -0.94997, -0.404481, -1.350041, 1.360461)
y=c(-1.89666, -2.27684, -0.459432, -0.498989, 3.005143, -3.375039)
lm.out<-lm(y~x)

newx = seq(min(x),max(x),by = 0.05)
conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence", level = 0.95)
plot(x, y, xlab="Binding affinity difference", ylab="Log2(Alt FPKM/REF FPKM)", main="Regression")
abline(lm.out, col="black")
matlines(newx, conf_interval[,2:3], col = "lightblue", lty=2)

#### FUS
x = c(-1.639152, 1.639152, -1.639152, 1.639152, -1.639152, -1.639152, -1.639152, -1.639152, 1.639152)
y = c(-1.517607, 2.342392, -1.031027, 0.703144, 0.403356, -1.186413, -0.678072, -1.186413, 1.234465)
lm.out<-lm(y~x)
newx = seq(min(x),max(x),by = 0.05)
conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence", level = 0.95)
plot(x, y, xlab="Binding affinity difference", ylab="Log2(Alt FPKM/REF FPKM)", main="FUS")
abline(lm.out, col="black")
matlines(newx, conf_interval[,2:3], col = "lightblue", lty=2)
summary(lm.out)
