########################### linear regression (motifbreak ) #############
### Fig5E
### RBPs for which changes in binding affinities were significantlly correlated with ASm6A imbalance
# setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/motifbreakR/by_ind/plot/")
#### 
name <- "TAF15"
x=c(1.639152, -1.639152, -1.639152, -1.639152, 1.639152, -1.639152)
y=c(-0.680382,1.234206,1.43758,1.974909,-0.750641,1.38822)
lm.out<-lm(y~x)
newx = seq(min(x),max(x),by = 0.05)
conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence", level = 0.95)
plot(x, y, xlab="Binding affinity difference", ylab="Log2(Alt FPKM/REF FPKM)", main=name)
abline(lm.out, col="black")
matlines(newx, conf_interval[,2:3], col = "lightblue", lty=2)

####
name <- "SRSF9"
x=c(-1.350041,0.943519,0.442079,1.639152)
y=c(0.825321,-0.219169,0.194945,-0.354301)
lm.out<-lm(y~x)
newx = seq(min(x),max(x),by = 0.05)
conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence", level = 0.95)
plot(x, y, xlab="Binding affinity difference", ylab="Log2(Alt FPKM/REF FPKM)", main=name)
abline(lm.out, col="black")
matlines(newx, conf_interval[,2:3], col = "lightblue", lty=2)
#### FUS
x = c(1.639152, 1.639152, -1.639152, -1.639152, -1.639152)
y = c(0.489204, 0.316284, -0.392753, -0.698876, -1.260652)
lm.out<-lm(y~x)
newx = seq(min(x),max(x),by = 0.05)
conf_interval <- predict(lm.out, newdata=data.frame(x=newx), interval="confidence", level = 0.95)
plot(x, y, xlab="Binding affinity difference", ylab="Log2(Alt FPKM/REF FPKM)", main="FUS")
abline(lm.out, col="black")
matlines(newx, conf_interval[,2:3], col = "lightblue", lty=2)
summary(lm.out)
#############


