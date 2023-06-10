######################## exomePeak ##########################
library(exomePeak)
setwd("/Data_Resource/m6A-seq/com_fetus_vs_adult")
#####
ip1 = "fetal_data/heart/ip/heart_IP_1-uniq.bam" 
ip2 = "fetal_data/heart/ip/heart_IP_2-uniq.bam"
ip3 = "fetal_data/heart/ip/heart_IP_3-uniq.bam"
input1 = "fetal_data/heart/input/heart_1_L6-uniq.bam"
input2 = "fetal_data/heart/input/heart_2_L7-uniq.bam"
input3 = "fetal_data/heart/input/heart_3_L6-uniq.bam"
IP_BAM=c(ip1, ip2, ip3)
INPUT_BAM=c(input1, input2, input3)
t_ip1 = "adult_data/heart/ip/CRR042290-uniq.bam"
t_ip2 = "adult_data/heart/ip/CRR055527-uniq.bam"
t_input1 = "adult_data/heart/input/CRR042291-uniq.bam"
t_input2 = "adult_data/heart/input/CRR055528-uniq.bam"
TREATED_IP_BAM=c(t_ip1, t_ip2)
TREATED_INPUT_BAM=c(t_input1, t_input2)
# comparison
exomepeak(GENOME="hg19", IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM, TREATED_IP_BAM=TREATED_IP_BAM, TREATED_INPUT_BAM=TREATED_INPUT_BAM, EXPERIMENT_NAME="heart")


########################################################### MeTDiff #############################################################
setwd("/data6/m6AQTL/data/diff_methylation/heart")

library(MeTDiff)

# in the real case, change the gtf to what you need
gtf <- "../../hg19_only_22/hg19.gtf"

ip1 <- 'fetus1.m6A.bam'
ip2 <- 'fetus2.m6A.bam'
ip3 <- 'fetus3.m6A.bam'
input1 <- 'fetus1.input.bam'
input2 <- 'fetus2.input.bam'
input3 <- 'fetus3.input.bam'
treated_ip1 <- 'adult1.m6A.bam'
treated_ip2 <- 'adult2.m6A.bam'
treated_input1 <- 'adult1.input.bam'
treated_input2 <- 'adult2.input.bam'
data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAABCklEQVR42p2ULQ7CMBiGISQIDAqHxKLxOwAX2Am4ACEkKAQ3IFkgbbd1XVJH0Hg0FskFsBhKi4L2XdJvb1K150nWn/frrDdbw3MVLFGoZecnsVwQY0y3QX4LUaVUDkZr3W+QX6xQCZWDyTI9bJCfnMsplYM5SDlGssjVw32jcjDHspxgub4zpkdUDoYV9azhGK5Zdh5QORgh6jmSWa5O9rJ7VA7vRFQL/PbrfRsuiIVWSHR/3IbzjqhKXYmQ5EpH5bxLVokrD5Jc2ajc/7HYsnxLAyRXMioXFM2VBUl+0WK4YFRY6IYkf1TEcMGws9AFSf6wi+GCcc0LJZHkj+sYDr31HdyufdttOD8f3G/pQdD2FdIAAAAASUVORK5CYII=
IP_BAM <- c(ip1,ip2,ip3)
INPUT_BAM <- c(input1,input2,input3)
TREATED_IP_BAM <- c(treated_ip1, treated_ip2)
TREATED_INPUT_BAM <- c(treated_input1, treated_input2)

metdiff(GENOME = "hg19", IP_BAM = IP_BAM,INPUT_BAM = INPUT_BAM, #GENE_ANNO_GTF=gtf,
        TREATED_IP_BAM = TREATED_IP_BAM,TREATED_INPUT_BAM=TREATED_INPUT_BAM,
        EXPERIMENT_NAME="heart")



##################################
############ RADAR ##############
library(RADAR)
# setwd("/Data_Resource/m6A-seq/com_fetus_vs_adult/RADAR_result/")

radar <- countReads(
  samplenames = c("stomachfetus1","stomachfetus2","stomachadult1","stomachadult2"), #1
  gtf = "/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_methylation/RADAR_result/hg19.gtf",
  bamFolder = "/Data_Resource/m6A-seq/com_fetus_vs_adult/bams/all_bams/", #2
  outputDir= "/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_methylation/RADAR_result/stomach",
  modification = "m6A",
  binSize = 50,
  strandToKeep = "both",  # Default is “opposite” for reverse-stranded library;set to “both” to ignore strand when counting the reads
  paired = FALSE,
  threads = 10,
  saveOutput = T
)
summary(radar)

radar <- normalizeLibrary( radar )
sizeFactors(radar)

radar <- adjustExprLevel( radar )
variable(radar) <- data.frame( group = c("Ctl","Ctl","Treated","Treated") )
radar <- filterBins( radar ,minCountsCutOff = 15)
radar <- diffIP_parallel(radar, thread = 10, fdrBy = "fdr")
radar <- reportResult( radar, cutoff = 0.05, Beta_cutoff = 0.5,threads = 10)
res <- results( radar )
write.table(res, "/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_methylation/RADAR_result/stomach/stomach_result.txt", sep="\t", quote = FALSE)


###################################### read RDS file ###############################
radar <- readRDS("/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_methylation/RADAR_result/brain/MeRIP_readCounts.RDS")
summary(radar)
radar <- normalizeLibrary(radar)
sizeFactors(radar)
radar <- adjustExprLevel( radar )
# brain
variable(radar) <- data.frame( group = c("Ctl","Ctl","Ctl","Treated","Treated","Treated","Treated","Treated","Treated") )
# heart
# variable(radar) <- data.frame( group = c("Ctl","Ctl","Ctl","Treated","Treated") )

radar <- filterBins( radar ,minCountsCutOff = 15)

# pca plot
top_bins <- extractIP(radar,filtered = T)[order(rowMeans( extractIP(radar,filtered = T) ),decreasing = T)[1:1000],]
plotPCAfromMatrix(top_bins,group = unlist(variable(radar)) )
#
radar <- diffIP_parallel(radar, thread = 30, fdrBy = "fdr", exclude = c("brainadult2","brainadult5","brainadult6")) # ~7min
radar <- reportResult( radar, cutoff = 0.001, Beta_cutoff = 0.7,threads = 30) # ~21min
# heatmap plot
plotHeatMap(radar)
#
res <- results( radar )
write.table(res, "/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_methylation/RADAR_result/brain/brain_result2.txt", sep="\t", quote = FALSE)


####################################### For all samples ############################
radar <- countReads(
  samplenames = c("brainadult1", "brainadult2", "brainadult3", "brainadult4", "brainadult5", "brainadult6", "brainfetus1", "brainfetus2", "brainfetus3", "heartadult1", "heartadult2", "heartfetus1", "heartfetus2", "heartfetus3", "liveradult1", "liverfetus1", "liverfetus2", "liverfetus3", "lungadult1", "lungadult2", "lungadult3", "lungadult4", "lungfetus1", "lungfetus2", "stomachadult1", "stomachadult2", "stomachfetus1", "stomachfetus2"), #1
  gtf = "/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_methylation/RADAR_result/hg19.gtf",
  bamFolder = "/Data_Resource/m6A-seq/com_fetus_vs_adult/bams/all_bams/", #2
  outputDir= "/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_methylation/RADAR_result/all_samples",
  modification = "m6A",
  binSize = 100,
  strandToKeep = "both", 
  paired = FALSE,
  threads = 10,
  saveOutput = T
)

summary(radar)
radar <- normalizeLibrary( radar )
sizeFactors(radar)
radar <- adjustExprLevel( radar )
variable(radar) <- data.frame( group = c("Treated", "Treated", "Treated", "Treated", "Treated", "Treated", "Ctl", "Ctl", "Ctl", "Treated", "Treated", "Ctl", "Ctl", "Ctl", "Treated", "Ctl", "Ctl", "Ctl", "Treated", "Treated", "Treated", "Treated", "Ctl", "Ctl", "Treated", "Treated", "Ctl", "Ctl") )
radar <- filterBins( radar ,minCountsCutOff = 15)

# pca plot
top_bins <- extractIP(radar,filtered = T)[order(rowMeans( extractIP(radar,filtered = T) ),decreasing = T)[1:1000],]
plotPCAfromMatrix(top_bins,group = unlist(variable(radar)) )
#
radar <- diffIP_parallel(radar, thread = 20, fdrBy = "fdr") # ~7min
# radar <- reportResult( radar, cutoff = 0.001, Beta_cutoff = 0.7,threads = 20) # ~21min
# # heatmap plot
# plotHeatMap(radar)
# 
# res <- results( radar )
# write.table(res, "/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_methylation/RADAR_result/all_samples/result.txt", sep="\t", quote = FALSE)
