library(ChIPseeker)
setwd("/home/galaxy/project/alleleSpecific_analysis/results/")
#m6a <- "ASE_analysis/readCount/Fisher_test/sig/bed/extend_bed/2kb/CRR055543.bed"
#con = "ASM/corr_m6A/corr_results/plot_profile/m6A_region/ASM_control.bed2"
#amr = "ASM/corr_m6A/corr_results/plot_profile/m6A_region/ASM.bed2"
###############
sig="ASE_analysis/readCount/Fisher_test/sig01/bed/extended_10kb/Spleen-3-2.bed"
#unsig="ASM/corr_m6A/corr_results/plot_profile/AMR_region/extend_bed/CRR055551_control.bed"
amr="ASM/allelicEpigenome/by_tissue/bed/extend_250bp/Spleen.bed"
# 
amr_peak <- readPeakFile(amr)
#tagMatrixList <- lapply(c(sig, unsig), getTagMatrix, windows=amr_peak)
#names(tagMatrixList) <- c("Allele-specific m6A", "Control")
m6a_peak <- readPeakFile(sig)
# tagMatrixList <- lapply(c(amr, con), getTagMatrix, windows=m6a_peak)
# names(tagMatrixList) <- c("AMR", "Control")
amrMatrix <- getTagMatrix(amr_peak, windows = m6a_peak)
#conMatrix <- getTagMatrix(con_peak, windows = m6a_peak)
#plotAvgProf(tagMatrixList, xlim=c(-1000, 1000),xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency")
plotAvgProf(amrMatrix, xlim=c(-10000, 10000),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

