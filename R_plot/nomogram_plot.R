library(pROC)
y_true = c(0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0)
y_predict = c(0.6060616554576405,0.16169185699231362,0.9971093143280013,0.8034536466575728,0.07374722987667176,0.3413087198383404,0.18392731200992293,0.1807230365644545,0.748678337579735,0.254130412385619,0.14539254269667673,0.2514809796650626,0.2289191206796813,0.17885915584767556,0.5085488127165505,0.058153306804457974,0.07054741462243307,0.17087561708957047,0.9655306101680667,0.07194206310965236,0.03377967866691865,0.4417878211170693,0.03722473239398826,0.0631916144987044,0.6082205875468354,0.7615779259833511,0.21793290031707552,0.2369139661561988,0.19018527571999205,0.22884311928888712,0.9315726150940291,0.2961111034049563,0.20496029294195683,0.06157951560097004)
roc = roc(y_true, y_predict)
roc
ci(roc)



library(rms)
setwd("/home/galaxy/project/logistic_regression/new")
in_file = "/home/galaxy/project/logistic_regression/new/nomogram_total.txt" # 
df = read.table(in_file, sep = "\t", header =TRUE)
head(df)
# colnames(df) <- c("c4123", "c428112", "c428422", "c428522", "Cesarean", "label")
y = df$Label
Age <- df$Age
Cesarean <- df$Cesarean
Abortion <- df$abortion
Ultrasound <- df$ultrasound.PA.
circScore <- df$circScore
#
Age <- factor(Age)
Ultrasound <- factor(Ultrasound)
Abortion <- factor(Abortion)
Cesarean <- factor(Cesarean)
circScore <- factor(circScore)
# BMI <- df$BMI
ddist <- datadist(Age, Cesarean, Abortion, Ultrasound, circScore)
options(datadist='ddist')
f <- lrm(y ~ Age+Cesarean+Abortion+Ultrasound+circScore)
nom <- nomogram(f, fun=plogis,
                fun.at=c(.001, .01, seq(.1,.9, by=.4), .99, .999),
                lp=F, funlabel="Risk")
plot(nom)
