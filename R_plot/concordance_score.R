data("faithful")
head(faithful)
## display the empirical distribution of the "wating"
y = faithful$waiting
hist(y, xlab = "waiting time (mn)")
## we clearly see 2 modes: the waiting time seem to be distributed either around 50mn or around 80mn



##############################################################
##############################################################
data("diabetes")
summary(diabetes)

