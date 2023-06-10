library(Sushi)

##################################################################################################
#   ip vs ip
#################################################################################################
adult_file <- "/Data_Resource/m6A-seq/com_fetus_vs_adult/bams/adult/ip/CRR055527-uniq.bedGraph"
adult <- read.table(adult_file, sep = "\t", header = FALSE)
head(adult)
adult$V4 <- log2(adult$V4)

fetal_file <- "/Data_Resource/m6A-seq/com_fetus_vs_adult/bams/fetal/ip/heart_IP_1-uniq.bedGraph"
fetal <- read.table(fetal_file, sep = "\t", header = FALSE)
head(fetal)
fetal$V4 <- log2(fetal$V4)

chrom <- "chr20"
chromstart <- 1
chromend <- 63025520
par(mfrow=c(2,1),mar=c(1,4,1,1))

plotBedgraph(adult,chrom,chromstart,chromend,transparency=.50,color=SushiColors(2)(2)[1])
axis(side=2,las=2,tcl=.2)
mtext("Adult Read Depth",side=2,line=1.75,cex=1,font=2)
labelgenome(chrom,chromstart,chromend,n=5,scale="Mb")

plotBedgraph(fetal,chrom,chromstart,chromend,
             transparency=.50,color=SushiColors(2)(2)[2],flip=TRUE,
             rescaleoverlay=TRUE)
axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),
     labels=-1*pretty(par("yaxp")[c(1,2)]))
mtext("Fetus Read Depth",side=2,line=1.75,cex=1,font=2)
labelgenome(chrom,chromstart,chromend,n=5,scale="Mb")

##############################################################################################
# ip vs input
##############################################################################################
ip_file <- "/Data_Resource/m6A-seq/com_fetus_vs_adult/bams/fetal/ip/heart_IP_1-uniq.bedGraph"
ip <- read.table(ip_file, sep = "\t", header = FALSE)
head(ip)
ip$V4 <- log2(ip$V4)

input_file <- "/Data_Resource/m6A-seq/com_fetus_vs_adult/bams/fetal/input/heart_1_L6-uniq.bedGraph"
input <- read.table(input_file, sep = "\t", header = FALSE)
head(input)
input$V4 <- log2(input$V4)

chrom <- "chr20"
chromstart <- 1
chromend <- 630255
par(mfrow=c(2,1),mar=c(1,4,1,1))

plotBedgraph(ip,chrom,chromstart,chromend,transparency=.50,color=SushiColors(2)(2)[1])
axis(side=2,las=2,tcl=.2)
mtext("IP Read Depth",side=2,line=1.75,cex=1,font=2)
labelgenome(chrom,chromstart,chromend,n=5,scale="Mb")

plotBedgraph(input,chrom,chromstart,chromend,
             transparency=.50,color=SushiColors(2)(2)[2],flip=TRUE,
             rescaleoverlay=TRUE)
axis(side=2,las=2,tcl=.2,at=pretty(par("yaxp")[c(1,2)]),
     labels=-1*pretty(par("yaxp")[c(1,2)]))
mtext("Input Read Depth",side=2,line=1.75,cex=1,font=2)
labelgenome(chrom,chromstart,chromend,n=5,scale="Mb")

