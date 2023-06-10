library("peer")
library("MatrixEQTL")

# library(devtools)
# install_github("InouyeLab/BootstrapQTL")
# 
# # CRAN Install:
# install.packages("BootstrapQTL")
setwd("/home/galaxy/project/m6AQTL/2019_10_10/transcriptome_snp/SNPiR/shell/all_samples/results/MatrixEQTL_result/")
#############################################################################################
######### Find hidden batch effects and other confounders from the expression data ##########
#############################################################################################
##################### quantile normalization #########################
df = read.table("bedtools_quantification.matrix")
head(df)
df_norm = normalize.quantiles(data.matrix(df))
row.names(df_norm) = row.names(df)
colnames(df_norm) = colnames(df)
head(df_norm)
write.csv(df_norm, "m6A_norm.matrix", quote = FALSE)
#######################################################################
# df_t.to_csv("bedtools_quantification_T.matrix", index=False, header=False)
expr = read.csv("m6A_norm.matrix")
dim(expr)
## including covariates which contribute to variability in the data ###
# covs = read.csv("Covariates_T.csv", header = FALSE)
# PEER_setCovariates(model, as.matrix(covs))
########
model = PEER()
PEER_setPhenoMean(model,as.matrix(expr))
# PEER_setAdd_mean(model, TRUE) # PEER recommends include an additional factor to account for the mean expression.
dim(PEER_getPhenoMean(model))

PEER_setNk(model,5) # NULL response means no error here
PEER_getNk(model)

PEER_update(model)
factors = PEER_getX(model)
factors[,2]
dim(factors)
factors <- t(factors)
row.names(factors) <- sprintf("peer_%d", c(seq(1:5)))
write.csv(factors, "peer_cov.csv", col.names = F, quote = F)
#######                  ##########
### other matrics about factors ### 
###################################
#weights = PEER_getW(model)
#dim(weights)
#precision = PEER_getAlpha(model)
#residuals = PEER_getResiduals(model)
#dim(residuals)
#plot(precision)

# Run the BootstrapQTL analysis
# eGenes <- BootstrapQTL(snps, gene, snpspos, genepos,
#                        n_bootstraps=200, n_cores=2,
#                        eGene_detection_file_name = "cis_eQTL_associations.txt",
#                        bootstrap_file_directory = "bootstrap_analyses/")


##############################################################################
############################ cis and trans QTL ###############################
##############################################################################
setwd("/home/galaxy/project/m6AQTL/2019_10_10/transcriptome_snp/SNPiR/shell/all_samples/results/MatrixEQTL_result")
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "SNP_rename.csv"
snps_location_file_name = "snpsloc.txt"

# Gene expression file name
expression_file_name = "m6A_norm.matrix"
gene_location_file_name = "m6Aloc.csv"

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "peer_cov.csv"

# Output file name
output_file_name_cis = "m6AQTL_cis.txt";
output_file_name_tra = "m6AQTL_trans.txt";

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.05; # 1e-5;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");
# Distance for local gene-SNP pairs
cisDist = 1e6;
## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = ",";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 3000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load m6A methylation data
gene = SlicedData$new();
gene$fileDelimiter = ",";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 3000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 0;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}
dim(cvrt)

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.csv(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot", # TRUE
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

# unlink(output_file_name_tra);
# unlink(output_file_name_cis);

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)
## Plot the Q-Q plot of local and distant p-values
plot(me)
####################
############
###############
########
####################
#####
model = PEER()
PEER_setNk(model,k)
PEER_setPhenoMean(model,as.matrix(expr_cov)) # Give PEER your expression data
PEER_setAdd_mean(model, TRUE)
PEER_setCovariates(model, as.matrix(covar4))  # PEER ask NxC matrix, where N=samples and C=covariates
PEER_setNmax_iterations(model, 100000)
PEER_update(model)
residuals = PEER_getResiduals(model)  # convert to GxN
colnames(residuals) = colnames(expr_cov)
rownames(residuals) = rownames(expr_cov)
# add mean
residuals = residuals + apply(expr_cov, 2, mean)
# convert to SliceData format
genes = SlicedData$new();
genes$CreateFromMatrix(residuals); # using your residuals as expression input
rm(residuals, model)


