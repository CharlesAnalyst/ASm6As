# https://rdrr.io/github/Jiaxin-Fan/ASEP/f/vignettes/introduction.Rmd
library(ASEP)
library(knitr)
library(DT)

setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASEP/")

dat = read.table("brain_adult_int.txt", header=T)
dat = dat[order(dat$id,dat$snp),]
dat = dat[,c("gene","id","snp","ref","total")]

ASE_detection(dat_all=dat,phased=FALSE,adaptive=TRUE,n_resample=10^3,parallel=TRUE,save_out=TRUE,name_out="brain_adult_results.txt")

dat$ref_condition = 'M1'
set.seed(6)