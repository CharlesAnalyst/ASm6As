library(Guitar)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
# import different data formats into a named list object.
# setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/Guita_plot/")
setwd("/Charles/project/ASm6A/ASm6A/guitar_plot/")
##########
# brain <- "Brainstem-5-3.bed" # "brain_IP_1.bed"
# heart <- "Heart-1-1.bed" #"heart_IP_1.bed"
# liver <- "Liver-4-2.bed" #"liver_IP_1.bed"
# lung <- "Lung-2-1.bed" #"lung_IP_4.bed"
# stomach <- "Stomach-4-2.bed"#"stomach_IP_4.bed"
# m6a_hg19 <- list(brain,heart,liver,lung,stomach)
# names(m6a_hg19) <- c("brain","heart","liver","lung","stomach")
# Build Guitar Coordinates
# txdb_file <- makeTxDbFromGFF("../../data/hg19_genome/hg19.gtf")
# txdb <- makeGuitarTxdb(txdb_file)
# GuitarPlot(stBedFiles=m6a_hg19, txTxdb = txdb)  # miscOutFilePrefix = "example"
# GuitarPlot(stBedFiles=m6a_hg19, txGenomeVer = "hg19", pltTxType = c("mrna"), enableCI = FALSE)
#headOrtail = TRUE,
#enableCI = FALSE,
#mapFilterTranscript = TRUE,
#pltTxType = c("mrna"),
#stGroupName = c("heart")
####################################################################
# asm6a <- "total.bed"
# xia <- "xia_hg19_peak.bed"
# yi <- "yi_hg19_peak.bed"
# zhang <- "zr_hg19_peak.bed"
# m6a_hg19 <- list(asm6a, xia, yi, zhang) # asm6a, xia, yi, zhang
# names(m6a_hg19) <- c("ASm6As", "xia","yi","zr")
#asm6a <- "ASm6A.bed"
m6a <- "total.bed"
m6a_hg19 <- list(m6a) # asm6a, xia, yi, zhang
names(m6a_hg19) <- c("ASm6As") # "ASm6As", "Fetal peak (Xiao)", "Adult peak (Liu)", "Adult peak (Zhang)"
# Build Guitar Coordinates
# txdb_file <- makeTxDbFromGFF("/home/galaxy/project/alleleSpecific_analysis/data/hg19_genome/hg19.gtf")
# txdb <- makeGuitarTxdb(txdb_file)
# GuitarPlot(stBedFiles=m6a_hg19, txdb = txdb, pltTxType = c("mrna"))  # miscOutFilePrefix = "example"
GuitarPlot(stBedFiles=m6a_hg19, txGenomeVer = "hg19", pltTxType = c("mrna"), enableCI = FALSE)
########################### for test ################################
# stBedFiles <- list(system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed",
#                                package="Guitar"))
# txdb_file <- system.file("extdata", "mm10_toy.sqlite",
#                          package="Guitar")
# txdb <- loadDb(txdb_file)
# guitarTxdb <- makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE)
# GuitarPlot(txTxdb = txdb,
#            stBedFiles = m6a_hg19,
#            miscOutFilePrefix = "example")
################
library(MeRIPtools)
# Plot Meta Gene of the peak you called.
#Here I use R package Guitar to plot the meta gene distribution of peaks. Guitar normalized the density of peaks on each feature (5UTR, CDS, 3UTR) by the length of the feature, which correlates to the expected occurence of peak by chance.
#To get start, we need to have the peak in bed12 format read into R. In previous peak calling step, if you used MeRIPtools to call joint peak:
Joint_peak <- reportJointPeak(TestRun, threads = 6)
#You already have your peak in variable Joint_peak. You can proceed to next step
plotMetaGene(Joint_peak,gtf = "~/Database/genome/hg38/hg38_UCSC.gtf")
# If you used other tools to call peak, you can first read peak file into R and then plot the meta gene:
peak <- read.table("path/to/peak.bed",sep = "\t", stringsAsFactors = F)
plotMetaGene(peak,gtf = "~/Database/genome/hg38/hg38_UCSC.gtf")
