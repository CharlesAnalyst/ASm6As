library("limma")
library("edgeR")
library("statmod")
library("DESeq2")

setwd("/Data_Resource/m6A-seq/com_fetus_vs_adult/m6A_abundance/quantification/norm_by_scaleFactor")


limma_voom_DE <- function(infile,a,b,method){
  infile = "bedtools_quantification.matrix"
  data = read.table(infile, header = T,row.names = 1) # Reading use input filename
  head(data)
  #data_tissue <- data[,1:9]  #brain
  #data_tissue <- data[,19:24]#lung
  # data_tissue <- data[,25:28]#stomach
  # data_tissue <- data[,15:18]#liver
  data_tissue <- data[,10:14]#heart
  # log_data <- log2(data_tissue+0.5)
  log_data <- data_tissue
  #head(log_data)
  #
  # brain
  # data_reorder <- data.frame(log_data$brainfetus1, log_data$brainfetus2, log_data$brainfetus3, log_data$brainadult1, log_data$brainadult2, log_data$brainadult3, log_data$brainadult4, log_data$brainadult5, log_data$brainadult6)
  # heart
  data_reorder <- data.frame(log_data$heartfetus1, log_data$heartfetus2, log_data$heartfetus3, log_data$heartadult1, log_data$heartadult2)
  # liver
  # data_reorder <- data.frame(log_data$liverfetus1, log_data$liverfetus2, log_data$liverfetus3, log_data$liveradult1)
  # lung
  #data_reorder <- data.frame(log_data$lungfetus1, log_data$lungfetus2, log_data$lungadult1, log_data$lungadult2, log_data$lungadult3, log_data$lungadult4)
  # stomach
  #data_reorder <- data.frame(log_data$stomachfetus1, log_data$stomachfetus2, log_data$stomachadult1, log_data$stomachadult2)
  rownames(data_reorder) <- rownames(log_data)
  head(data_reorder)
  group = as.factor(c(rep("fetus",3), rep("adult",2))) # Giving number of replicate information for each condition.
  group
  # nf <- calcNormFactors(data,method = "TMM") # TMM Normalization
  design <- model.matrix(~group) # Base on conditions and plicate information designing the info for your count matrix
  y <- voom(data_reorder, design, plot=TRUE) # , lib.size=colSums(data_heart)*nf
  fit <- lmFit(y,design) # Fitting to linear model for voom output
  fit <- eBayes(fit)
  # toptable = topTable(fit,coef=ncol(design),number=20000,adjust.method="BH") # Getting toptable significant hits
  toptable = topTable(fit,coef=ncol(design),number=nrow(log_data),adjust.method="fdr",sort.by="B",resort.by=NULL,p.value=0.05,lfc=0)
  # Results Writing to file
  write.table(toptable,"/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_methylation/limmaVoom_result/heart_diff2.txt",sep = "\t", quote = FALSE)
}