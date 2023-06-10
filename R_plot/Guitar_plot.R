library(Guitar)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
# import different data formats into a named list object.
setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/Guita_plot/")
##########
brain <- "Brainstem-5-3.bed" # "brain_IP_1.bed"
heart <- "Heart-1-1.bed" #"heart_IP_1.bed"
liver <- "Liver-4-2.bed" #"liver_IP_1.bed"
lung <- "Lung-2-1.bed" #"lung_IP_4.bed"
stomach <- "Stomach-4-2.bed"#"stomach_IP_4.bed"
m6a_hg19 <- list(brain,heart,liver,lung,stomach)
names(m6a_hg19) <- c("brain","heart","liver","lung","stomach")
# Build Guitar Coordinates
# txdb_file <- makeTxDbFromGFF("../../data/hg19_genome/hg19.gtf")
# txdb <- makeGuitarTxdb(txdb_file)
# GuitarPlot(stBedFiles=m6a_hg19, txTxdb = txdb)  # miscOutFilePrefix = "example"
GuitarPlot(stBedFiles=m6a_hg19, txGenomeVer = "hg19", pltTxType = c("mrna"), enableCI = FALSE)
#headOrtail = TRUE,
#enableCI = FALSE,
#mapFilterTranscript = TRUE,
#pltTxType = c("mrna"),
#stGroupName = c("heart")
####################################################################
##################### single file ##################################
asm6a <- "total.bed"
m6a_hg19 <- list(asm6a)
names(m6a_hg19) <- c("Common ASm6As")
# Build Guitar Coordinates
# txdb_file <- makeTxDbFromGFF("../../data/hg19_genome/hg19.gtf")
# txdb <- makeGuitarTxdb(txdb_file)
# GuitarPlot(stBedFiles=m6a_hg19, txTxdb = txdb)  # miscOutFilePrefix = "example"
GuitarPlot(stBedFiles=m6a_hg19, txGenomeVer = "hg19", pltTxType = c("mrna"), enableCI = FALSE)


########################### for test ################################
stBedFiles <- list(system.file("extdata", "m6A_mm10_exomePeak_1000peaks_bed12.bed",
                               package="Guitar"))
txdb_file <- system.file("extdata", "mm10_toy.sqlite",
                         package="Guitar")
txdb <- loadDb(txdb_file)
guitarTxdb <- makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE)
GuitarPlot(txTxdb = txdb,
           stBedFiles = m6a_hg19,
           miscOutFilePrefix = "example")
