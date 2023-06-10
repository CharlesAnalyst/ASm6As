library(ggplot2)
library("extrafont")
font_import()

#########################################################################
# in_file = "/Charles/project/ASm6A/plot/fig1c.txt"
# df <- read.table(in_file, sep = "\t", header = TRUE)
# print(head(df))
# df$Sample <- factor(df$Sample,levels = df$Sample,ordered = TRUE)
# color_list <- c("#ff6666","#ff9933","#53c653","#ffcc00","#6699ff", 'lightgrey')
# p <- ggplot(data=df, aes(x=Sample,y=Number, fill=Tissue)) + 
#   geom_bar(stat="identity", width=0.8) + 
#   coord_flip() +  
#   scale_fill_manual(values = color_list) +
#   xlab("") + 
#   ylab("ASm6A number") + 
#   theme_bw() + 
#   theme(
#     axis.title.x = element_text(size = 12, family = "Arial", color = "black"),
#     # axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
#     axis.ticks.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.text.x = element_text(family = "Arial", color = "black", size=12),
#     panel.border = element_blank(),
#     panel.grid = element_blank(),
#     axis.line = element_line(colour = "black")
#   )
# p
# ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig1c.pdf",
#        p, dpi=600, width = 3.27, height = 5.13)
##############################################################
##############################################################
## https://github.com/holtzy/R-graph-gallery/blob/master/299-circular-stacked-barplot.Rmd
## https://rpubs.com/Grady/630513
############################ circle barplot ##################
in_file = "/Charles/project/ASm6A/plot/fig1c.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
print(head(df))
df$Sample <- factor(df$Sample, levels = df$Sample, ordered = TRUE)
color_list <- c("#ff6666","#ff9933","#53c653","#ffcc00","#6699ff",'lightgrey')
#######
# This section prepare a dataframe for labels
# Get the name and the y position of each label
label_data <- df
label_data$id <- seq(1, nrow(df))
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360*(label_data$id-0.5)/number_of_bar  # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse(angle < -90, 1, 0)
# flip angle to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
#######
base_data <- data.frame(
  group = c("fetal","adult","fetal","adult","fetal","adult","fetal","adult","fetal","adult","fetal","adult"),
  start = c(1,4,16,19,24,27,30,33,37,39,46,48),
  end = c(3,15,18,23,26,29,32,36,38,45,47,50),
  angle <- c(79.2, 21.6, -28.8, -57.6, -86.4, 72.0, 50.4, 28.8, 0.0, -28.8, -57.6, -79.2) # extract angle from label_data (最中间那个angle)
)
base_data$title <- rowMeans(base_data[,c('start', 'end')], na.rm=TRUE)
head(base_data)

base_data2 <- data.frame(
  group = c("Brain", "Heart", "Kidney", "Liver", "Lung", "Muscle"),
  start = c(1,16,24,30,37,46),
  end = c(15,23,29,36,45,50),
  angle <- c(68, 16, -16, -44, -78, 10)
)
base_data2$title <- rowMeans(base_data2[,c('start', 'end')], na.rm=TRUE)
head(base_data2)
######
p <- ggplot(data=df, aes(x=Sample,y=Number,fill=Tissue)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = color_list) +
  ylim(-5000, 10000) + 
  theme_bw() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    # panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = unit(rep(-2,4), "cm")
  ) + 
  coord_polar(start = 0) + 
  # geom_text(data=label_data, aes(x=id, y=Number+26, label=Sample, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=3.5, angle= label_data$angle, inherit.aes=FALSE) +
# Add base line information
  geom_segment(data=base_data, aes(x=start, y=-150, xend=end, yend=-150),colour="black",alpha=0.8,size=0.6,inherit.aes=FALSE)  +
  geom_text(data=base_data,aes(x=title,y=-250,label=group),hjust=c(1,1,1,1,1,1,1,1,1,1,0,0),colour="black",alpha=0.8,size=4,fontface="bold", angle= base_data$angle, inherit.aes = FALSE) +
  geom_segment(data=base_data2, aes(x=start, y=-2250, xend=end, yend=-2250, colour=group),alpha=0.8,size=0.6,inherit.aes=FALSE)  + scale_colour_manual(values=color_list) + 
  geom_text(data=base_data2,aes(x=title,y=-2650,label=group),hjust=c(1,1,1,0,0,0),colour="black",alpha=0.8,size=4,fontface="bold", inherit.aes = FALSE) # angle= base_data2$angle, 


p
ggsave("/Charles/project/ASm6A/plot/fig1c.pdf",p, dpi=600, width = 9.86, height = 11)
#################################################################################################

##############################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/e_fig1a.txt"
df <- read.table(in_file, header = TRUE)
head(df)
p <- ggplot(df, aes(x=Proportion, color=Proportion)) +
  geom_histogram(fill="#009E73", bins = 50) + ylab("Frequency") + 
  scale_x_continuous(limits = c(0, 0.55),breaks=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6)) + theme_bw() + 
  theme(
    axis.title.x = element_text(size = 14, family = "Arial", color = "black"),
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
    axis.text.x = element_text(family = "Arial", color = "black", size=10),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"))
p
loadfonts() ############!!!!!!!!!!!!!!
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/e_fig1a.pdf", p, dpi=600, width = 15, height = 10.4, units = "cm")
###############################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/e_fig1c.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
print(head(df))

p <- ggplot(data=df, aes(x=Category, y=Proportion, fill=Category)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#999999", "#009E73"))+
  theme_bw() + labs(x="", y="Fraction") + guides(fill=FALSE) +
  theme(
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=10),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")
  ) + 
  scale_x_discrete(limits=c("ASm6A", "ASE"))
p
loadfonts()
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/e_fig1c.pdf", 
       p, dpi=600, width = 8, height = 8)
##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/VEF/count_results.txt"
df <- read.table(in_file)
head(df)
colnames(df) <- c("class", "n", "prop")
count.data <- df
# Add label position
# count.data <- count.data %>%
#   arrange(desc(class)) %>%
#   mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data
#
count.data$class <- factor(count.data$class, levels = count.data$class)
p <- ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  # geom_text(aes(y = lab.ypos, label = prop), color = "black")+
  scale_fill_brewer(palette="Set3") +
  theme_void()+
  theme(legend.text=element_text(size=14))
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/VEF/pie2.pdf", p, dpi=600, width = 15, height = 10.4, units = "cm")
##########################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig2a.txt"
df <- read.table(in_file, sep = "\t")
print(df)
df$V1 = factor(df$V1, levels = df$V1)

p <- ggplot(data=df, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.8)+
  scale_fill_manual(values=c("#D55E00", "#D55E00", "#D55E00", "#009E73", "#009E73"))+
  theme_bw() + guides(fill=FALSE)+ labs(x="", y="Log2 OR") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 60, hjust = 1, family = "Arial", color = "black", size=10),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig2a.pdf", 
       p, dpi=600, width = 3.6, height = 4.13)
#############################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/e_fig2a.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
head(df)
# Add label position
df <- df %>%
  arrange(desc(Category)) %>%
  mutate(lab.ypos = cumsum(Fraction) - 0.5*Fraction)
#
df$Category <- factor(df$Category, levels = df$Category)
p <- ggplot(df, aes(x = "", y = Fraction, fill = Category)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  #geom_text(aes(y = lab.ypos, label = Fraction), color = "black")+
  scale_fill_brewer(palette="Set3") +
  theme_void()+
  theme(legend.text=element_text(size=14))
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/e_fig2a.pdf", p, dpi=600, width = 15, height = 10.4, units = "cm")
#######################################################################
library(ggpubr)
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig2b.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
df$OR <- log2(df$OR)
df <- df[order(df$OR),]
df$RBP = factor(df$RBP, levels = df$RBP)
print(head(df))

df$Category<-factor(df$Category,levels=c('Writer','Reader','Regulator','Cofactor','Unknown'))
p <- ggbarplot(df, x = "RBP", y = "OR",
               fill = "Category",
               color = "White",                               
               palette = c("#D55E00", "#E69F00", "#00AFBB", "#47b8e0", "#999999"),
               # sorting = "descending",                       
               # add = "segments",                             
               xlab = "",
               rotate = TRUE,                                  
               group = "Category",
               ggtheme = theme_pubr()) +
  theme(legend.position = "right")
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig2b.pdf", 
       p, dpi=600, width = 7, height = 9)
###########################################################################################

enriched_rbps = c("YTHDF1", "UCHL5", "SND1", "YTHDF2", "ZNF622", "GRWD1", "EIF3H", "FXR2", "YBX3", "DDX55", "PPIG", "RPS3", "RBM15", "UPF1", "LARP4", "IGF2BP1", "EIF4A3", "FMR1", "METAP2", "TRA2A", "SUB1", "IGF2BP3", "YTHDF3", "IGF2BP2", "DDX24", "LSM11", "LIN28A", "ATXN2", "CAPRIN1", "NCBP3", "YTHDC1", "DDX6", "NOL12", "MOV10", "RBM15B", "PUM2", "DDX3X", "FXR1", "MBNL2", "SRRM4", "SRSF9", "RBM47", "EIF3D", "CPSF7", "TIAL1", "CPSF6", "LIN28B", "SERBP1", "BUD13", "RTCB", "NUDT21", "EIF3G", "GRSF1", "TIA1", "ZC3H7B", "MTPAP", "NOP58", "SRSF1", "ALKBH5", "GEMIN5", "YWHAG", "FUBP3", "FUS", "FAM120A", "EIF3B", "PCBP2", "XRCC6", "CPSF1", "FASTKD2", "SLTM", "AARS", "FIP1L1", "EIF3A", "TAF15", "ZRANB2", "HNRNPK", "DHX30", "FTO", "CPSF3", "AKAP8L", "SAFB2", "NCBP2", "TARDBP", "DGCR8", "AGGF1", "CELF2", "DROSHA", "EWSR1", "SF3B1", "CPSF4", "RBFOX2", "FBL", "TBRG4", "RBM22", "XPO5", "PTBP1", "SSB", "CSTF2T", "TNRC6A", "NONO", "TROVE2", "RPS5", "GTF2F1", "LARP7", "CSTF2", "SF3B4", "DDX59", "PRPF8", "HNRNPC", "NKRF", "WDR33", "CPSF2", "SRSF7", "U2AF2", "ELAVL1", "EFTUD2", "HNRNPH1", "HNRNPU", "XRN2", "ILF3", "SMNDC1")
cofactor_rbps = c("ABHD15", "ACAT1", "ACTB", "ACTC1", "ACTR2", "ACTR3", "AFMID", "AGPAT1", "AHNAK", "AKAP8L", "ALB", "ALDOA", "ALPI", "AMBRA1", "ANXA1", "ANXA2", "ANXA5", "ANXA6", "APEX1", "ARFIP2", "ARG1", "ARID5B", "ARNTL2", "ATAD3C", "ATG16L1", "ATP2A1", "ATP5A1", "ATP5B", "BAG1", "BAG2", "BAG3", "BAG5", "BAG6", "BANF1", "BCAP31", "BCL6", "BLMH", "BRIX1", "BRMS1", "BTN1A1", "BYSL", "C17orf59", "C18orf63", "C19orf68", "C1orf167", "C1orf68", "C1QBP", "C4A", "CA4", "CAD", "CALD1", "CALR", "CANX", "CAPZA2", "CASP14", "CAT", "CBLL1", "CBX3", "CCAR2", "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT7", "CCT8", "CD302", "CDK9", "CDKN2A", "CENPA", "CEP19", "CFL1", "CLTCL1", "CLU", "COL14A1", "CORO1C", "COX2", "COX5A", "CPS1", "CPSF6", "CS", "CSTA", "CTSD", "CUL2", "DCD", "DDB1", "DDB2", "DDI2", "DDOST", "DDX1", "DDX17", "DDX21", "DDX39A", "DDX39B", "DDX3X", "DDX3Y", "DDX46", "DDX5", "DEK", "DHX9", "DMBT1", "DNAJA1", "DNAJA2", "DNAJB1", "DNAJB11", "DNAJB4", "DNAJB6", "DNAJC7", "Dppa3", "DSC1", "DSG1", "DSP", "DYNLL1", "EBNA1BP2", "EEF1A1", "EEF1A2", "EEF1D", "EEF1G", "EEF2", "EIF2S2", "EIF3B", "EIF3G", "EIF3H", "EIF3I", "EIF3J", "EIF3M", "EIF4A1", "EIF4B", "EIF4E", "EIF5A", "ELAVL1", "EMD", "ENO1", "ERICH3", "ESR2", "ETFB", "EWSR1", "FAM172A", "FAM175A", "FASN", "FBXO38", "FLG", "FLG2", "FLNA", "FN1", "FSCN1", "FUS", "FZD8", "GANAB", "GAPDH", "GNB2L1", "GNL3", "GOT2", "GTF2I", "GTF3C4", "GTF3C5", "H2AFY", "H2AFZ", "H3F3C", "HDAC1", "HDAC2", "HIST1H1C", "HIST1H2AH", "HIST1H2BG", "HIST1H2BK", "HIST1H4A", "HIST2H3A", "HIST4H4", "HMGB1P1", "HNRNPA1", "HNRNPA2B1", "HNRNPA3", "HNRNPAB", "HNRNPC", "HNRNPCL1", "HNRNPD", "HNRNPDL", "HNRNPF", "HNRNPH1", "HNRNPK", "HNRNPM", "HNRNPR", "HNRNPU", "HPS6", "HRNR", "HSD17B10", "HSP90AA1", "HSP90AB1", "HSP90B1", "HSPA1A", "HSPA4", "HSPA4L", "HSPA5", "HSPA7", "HSPA8", "HSPA9", "HSPBP1", "HSPD1", "HSPE1", "HSPH1", "HTATSF1", "HUWE1", "IFI16", "IFT20", "IGHG1", "IGHG2", "IGKV1-5", "IKZF1", "ILF2", "ILF3", "IMPDH2", "IQGAP1", "JUP", "KCNH1", "KHDRBS3", "KHSRP", "KIAA1429", "KPRP", "KRT1", "KRT10", "KRT14", "KRT16", "KRT17", "KRT18", "KRT2", "KRT32", "KRT4", "KRT5", "KRT6A", "KRT6B", "KRT7", "KRT72", "KRT77", "KRT78", "KRT8", "KRT80", "KRT82", "KRT84", "KRT9", "KTN1", "KXD1", "LDHA", "LDHB", "LEFTY1", "LGALS7", "LIMA1", "LMNA", "LMNB1", "LRPPRC", "LRRC59", "LSM4", "LTF", "LYZ", "MAGEA11", "MATR3", "MCM3", "MED23", "METTL14", "METTL3", "MLF2", "MMS19", "MOV10", "MPP1", "MSH2", "MSH6", "MSN", "MTDH", "MYBBP1A", "MYC", "MYH9", "MYL6", "NACA", "NANOG", "NCBP1", "NCCRP1", "NCL", "NDUFAB1", "NFIB", "NGB", "NME1", "NOC2L", "NOLC1", "NONO", "NOP2", "NOP56", "NOP58", "NPM1", "NSUN2", "NUDT21", "NUP210", "OLFM2", "P4HB", "PAICS", "PARP1", "PCBP1", "PCBP2", "PDCD10", "PDHA1", "PDIA3", "PDIA6", "PFKP", "PHB", "PHB2", "PHGDH", "PIGB", "PLEC", "PLEKHA4", "PLG", "PMPCA", "PMPCB", "POLD1", "PPHLN1", "PPIA", "PPIB", "PPP1CC", "PRAME", "PRDX1", "PRDX2", "PRDX3", "PRKDC", "PRMT1", "PRMT5", "PRPF31", "PRPF8", "PRSS1", "PRSS3", "PSIP1", "PSMA1", "PSMC1", "PSMC3", "PSMC4", "PSMC5", "PSMC6", "PSMD13", "PTBP1", "PUF60", "PYHIN1", "RALY", "RAPGEF4", "RBBP7", "RBM14", "RBM15", "RBM15B", "RBMX", "RECQL4", "ROBO2", "RPL10", "RPL11", "RPL12", "RPL13", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL22", "RPL23", "RPL24", "RPL27", "RPL3", "RPL31", "RPL35", "RPL35A", "RPL36", "RPL38", "RPL4", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPLP0", "RPLP2", "RPN1", "RPN2", "RPS11", "RPS13", "RPS14", "RPS15A", "RPS16", "RPS18", "RPS2", "RPS25", "RPS26", "RPS27", "RPS3", "RPS3A", "RPS4X", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPSA", "RSL1D1", "RTN4", "RUVBL1", "RUVBL2", "RYR1", "SAMD9", "SDHA", "SEC22B", "SEL1L", "SERBP1", "SERPINB12", "SERPINH1", "SET", "SF3A3", "SF3B2", "SF3B3", "SF3B4", "SFPQ", "SH3BP5L", "SH3GL2", "SH3KBP1", "SHMT2", "SIRT6", "SLC25A3", "SLC25A5", "SLC25A6", "SLC3A2", "SLIRP", "SMAD2", "SMAD3", "SNF8", "SNRNP200", "SNRNP70", "SNRPB", "SNRPB2", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNW1", "SNX15", "Sod1", "SPATA5L1", "SPTAN1", "SPTBN1", "SRP14", "Srsf1", "SRSF1", "SRSF10", "SRSF2", "SRSF3", "SRSF6", "SRSF7", "SRSF9", "SSR1", "STIP1", "STT3A", "STUB1", "SUB1", "SYNCRIP", "TAB1", "TACC3", "TARDBP", "TBC1D10C", "TBL1XR1", "TBR1", "TCEB1", "TCEB2", "TCF12", "TCP1", "TFAM", "THAP11", "THRAP3", "TMEM200B", "TMPO", "TNIP1", "TNIP2", "TOP2A", "TP53BP1", "TPM1", "TPM3", "TPM4", "TRAP1", "TRAPPC6A", "TRIM25", "TRIM28", "TRIM4", "TSEN34", "TUBA1C", "TUBA4A", "TUBAL3", "TUBB", "TUBB2A", "TUBB4A", "TUBB4B", "TUFM", "TXN", "TXNDC5", "U2AF2", "U2SURP", "UBC", "UBR5", "UCKL1", "USP7", "VAPA", "VCP", "VDAC1", "VDAC2", "VIM", "VPS4A", "VPS52", "VSIG8", "WDR77", "WRNIP1", "WT1", "WTAP", "XRCC5", "XRCC6", "XRN2", "YBX1", "YWHAQ", "ZC3H13", "ZNF239", "ZNF572", "ZNF645", "ZRANB2", "ZSCAN25")
wjk_rbps = c("AKAP8L", "BUD13", "CAPRIN1", "CSTF2T", "DDX24", "DDX55", "ELAVL1", "GPKOW", "GTF2F1", "ILF3", "LSM11", "MOV10", "MSI1", "NCBP2", "RBM10", "RBM22", "RC3H1", "SF3B4", "SFPQ", "SLBP", "SUGP2", "TAF15", "TARDBP", "TRA2A", "UPF1", "YWHAG", "RBM4", "RBM41", "SNRNP70", "PABPC3", "RALY", "ZCRB")




########################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/e_fig2a.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
# df<- df[seq(dim(df)[1],1),]
print(df)
df$RBP = factor(df$RBP, levels = df$RBP)

p <- ggplot(data=df, aes(x=RBP, y=Number, fill=Enriched)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.8) +
  scale_fill_manual(values=c("#0072B2","#D55E00"))+
  theme_bw() + labs(x="", y="Number of ASm6As overlapping with ASB") + # + guides(fill=FALSE)
  theme(
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=10),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/e_fig2a.pdf", 
       p, dpi=600, width = 6.2, height = 2.57)

#####################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/e_fig2b.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
print(head(df))

p <- ggplot(data=df, aes(x=Category, y=Proportion, fill=Category)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#0072B2","#D55E00"))+
  theme_bw() + labs(x="", y="Fraction") + # + guides(fill=FALSE)
  theme(
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=10),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/e_fig2b.pdf", 
       p, dpi=600, width = 4.2, height = 6.8)

##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/e_fig2c.txt"
df <- read.table(in_file, sep = "\t")
print(head(df))

#
df$V1 <- factor(df$V1,levels = c('Motif loss allele','Motif gain allele'),ordered = TRUE)
p <- ggplot(data=df, aes(x=V1, y=V2, fill=V1)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#0072B2","#D55E00"))+
  theme_bw() + labs(x="", y="Log10(alt/ref)") + guides(fill=FALSE) +
  theme(
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=10),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
    # legend.title = element_blank()
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/e_fig2c.pdf", 
       p, dpi=600, width = 2.6, height = 4.5)

##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig3a.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
print(head(df))

#
df$Category <- factor(df$Category,levels = c('High-m6A allele','Low-m6A allele'),ordered = TRUE)
p <- ggplot(data=df, aes(x=Category, y=Allele.expression, fill=Category)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#0072B2","#D55E00"))+
  theme_bw() + labs(x="", y="Log10 FPKM") + guides(fill=FALSE) +
  theme(
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=12),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
    # legend.title = element_blank()
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig3a.pdf", 
       p, dpi=600, width = 5.2, height = 9.1)
###########################################################################
# file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig3b.txt"
file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig3b_2.txt"
df = read.table(file, header = TRUE)
head(df)

ggplot(df, aes(x = Proportion)) +
  stat_ecdf(aes(color = Category), geom = "step", size = 0.8) +
  scale_color_manual(values = c("orangered1", "steelblue2"))+
  labs(y = "Frequency") + theme_classic() +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig3d.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
print(head(df))

df$Region <- factor(df$Region,levels = df$Region,ordered = TRUE)
p <- ggplot(data=df, aes(x=Region, y=rho, group=1)) +
            geom_line(color="#0072B2") +
            geom_point() + ylim(0, -0.8) + 
            theme_bw() + labs(x="", y="Spearman correlation coefficient (rho)") +
  theme(
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=12),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")
    # legend.title = element_blank()
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig3d.pdf", 
       p, dpi=600, width = 5.36, height = 5)
##################################################################################
library(patchwork)

in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4ac.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
print(head(df))

p <- ggplot(df, aes(x = log10(Value), y = Tissue)) + 
  geom_point(aes(colour = Tissue, shape=Category), size=6) + 
  ylab("") + xlab("") + 
  guides() +
  geom_path(aes(group = Category), size=0.5, color="grey") + 
  scale_colour_manual(values = c("#ff6666","#ff9933","#53c653","#ffcc00","#6699ff")) + 
  theme_bw() + 
  theme(
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.x = element_text(hjust = 1, family = "Arial", color = "black", size=12),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig4ac.pdf", 
       p, dpi=600, width = 7.5, height = 6.5)
##################################################################################
# in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4a.txt"
# df <- read.table(in_file, sep = "\t") # , header = TRUE
# print(head(df))
# 
# df$V1 <- factor(df$V1,levels = df$V1,ordered = TRUE)
# color_list <- c("#ff6666","#ff9933","#53c653","#ffcc00","#6699ff")
# p <- ggplot(data=df, aes(x=V1, y=V2, fill=V1)) +
#   geom_bar(stat = "identity", width = 0.9) +
#   scale_fill_manual(values = color_list) +
#   theme_bw() + labs(x="", y="ASm6A number") + guides(fill=FALSE) +
#   theme(
#     axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
#     axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
#     axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=12),
#     panel.border = element_blank(),
#     panel.grid = element_blank(),
#     axis.line = element_line(colour = "black")
#     # aspect.ratio = 2/1
#   )
# p
# ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig4a.pdf", 
#        p, dpi=600, width = 5, height = 6)
##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4b.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
print(df)
df$Description <- factor(df$Description,levels = df$Description,ordered = TRUE)
color_list <- c("#ff6666","#ff9933","#53c653","#ffcc00","#6699ff")
p <- ggplot(data=df, aes(x=Description,y=X.log10p.adjust, fill=tissue)) + 
  geom_bar(stat="identity", width=0.8) + 
  coord_flip() +  
  scale_fill_manual(values = color_list) +
  xlab("") + 
  ylab("-log10p.adjust") + 
  theme_bw() + 
  theme(
        axis.title.x = element_text(size = 14, family = "Arial", color = "black"),
        axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
        axis.text.x = element_text(family = "Arial", color = "black", size=10),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black")
      )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig4b.pdf",
       p, dpi=600, width = 8, height = 10, units = "in")
##################################################################################
# in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4c.txt"
# df <- read.table(in_file, sep = "\t") # , header = TRUE
# print(head(df))
# 
# df$V1 <- factor(df$V1,levels = df$V1,ordered = TRUE)
# p <- ggplot(data=df, aes(x=V1, y=V2, fill=V1)) +
#   geom_bar(stat = "identity") +
#   scale_fill_brewer(palette = "Set2") +
#   theme_bw() + labs(x="", y="ASm6A number") + guides(fill=FALSE) +
#   theme(
#     axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
#     axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
#     axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=12),
#     panel.border = element_blank(),
#     axis.line = element_line(colour = "black")
#     # legend.title = element_blank()
#   )
# p
# ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig4c.pdf", 
#        p, dpi=600, width = 5.2, height = 6.2)
##################################################################################
library(ggpubr)
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4d.txt"
df <- read.table(in_file, sep = "\t") # , header = TRUE
colnames(df) <- c("stage","region","number")

df$number <- log2(df$number)
print(head(df))

p <- ggdotchart(df, x = "region", y = "number",
           color = "stage",                                # Color by groups
           # palette = c("#00AFBB", "#E7B800"), # Custom color palette
           # sorting = "descending",                       # Sort value in descending order
           # add = "segments",                             # Add segments from y = 0 to dots
           xlab = "",
           rotate = TRUE,                                  # Rotate vertically
           group = "stage",                                # Order by groups
           dot.size = 14,                                  # Large dot size
           label = round(df$number),                       # Add mpg values as dot labels
           font.label = list(color = "white", size = 22, 
                             vjust = 0.5),                 # Adjust label parameters
           ggtheme = theme_pubr()                          # ggplot2 theme
)
p
# df$V1 <- factor(df$V1,levels = df$V1,ordered = TRUE)
# p <- ggplot(data=df, aes(x=V2, y=V3, fill=V1)) +
#   geom_bar(stat = "identity", position=position_dodge()) +
#   scale_fill_brewer(palette = "Set2") +
#   theme_bw() + labs(x="", y="ASm6A number") + guides(fill=FALSE) +
#   theme(
#     axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
#     axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
#     axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=12),
#     panel.border = element_blank(),
#     axis.line = element_line(colour = "black")
#     # legend.title = element_blank()
#   ) + scale_x_discrete(limits=c("promoter-TSS", "5UTR", "exon", "intron", "3UTR", "TTS", "Intergenic", "non-coding"))
# p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig4d.pdf", 
       p, dpi=600, width = 8, height = 4)
##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/e_fig4a.txt"
df <- read.table(in_file, sep = "\t") # , header = TRUE
print(head(df))

# df$V1 <- factor(df$V1,levels = df$V1,ordered = TRUE)
p <- ggplot(data=df, aes(x=V2, y=V3, fill=V1)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() + labs(x="", y="ASm6A number") + guides(fill=FALSE) +
  theme(
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial", color = "black", size=12),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
    # legend.title = element_blank()
  ) + scale_x_discrete(limits=c("protein_coding", "antisense", "pseudogene", "processed_transcript", "sense_overlapping", "lincRNA", "3prime_overlapping_ncrna", "sense_intronic", "TR_C_gene"))
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/e_fig4a.pdf", 
       p, dpi=600, width = 8.3, height = 4.2)
########################################################################
library(ggpubr)

color_list <- c("#ff6666","#ff9933","#53c653","#ffcc00","#6699ff")
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig5a.txt"
df <- read.table(in_file, header = TRUE)
# df$Tissue <- factor(df$Tissue, levels = df$Tissue, ordered = FALSE) # 
head(df)
df<- df[seq(dim(df)[1],1),]

p <- ggdotchart(df, x = "Tissue", y = "OR",
                color = "Tissue",                                
                palette = color_list, 
                sorting = "descending",
                # sort.val = "OR",
                # add = "segments",                       
                rotate = TRUE,
                # group = "Group", 
                dot.size = 8,
                xlab = "",
                legend = "none",
                #label = round(df$OR),
                #font.label = list(color = "white", size = 14, vjust = 0.5), 
                ggtheme = theme_pubr())
p
# p2 <- ggpar(p, ylim = c(0.5, 2.3))
# p2
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig5a.pdf", p, dpi=600, width = 3.8, height = 5.8)
##########################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4b.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
print(df)
color_list <- c("#ff6666","#ff9933","#53c653","#ffcc00","#6699ff", 'grey')
p <- ggplot(data=df, aes(x=Tissue, y=OR, fill=Tissue)) +
  geom_bar(stat="identity", width=0.8) +
  coord_flip() +
  scale_fill_manual(values = color_list) +
  xlab("") +
  ylab("Odds ratio") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14, family = "Arial", color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(family = "Arial", color = "black", size=10),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black")
)
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig4b.pdf",
       p, dpi=600, width = 6, height = 4, units = "in")
##########################################################################
library(EnhancedVolcano)
file <- "/home/galaxy/project/alleleSpecific_analysis/results/ASE_analysis/readCount/Fisher_test/alleleRatio/Lung-4-2.txt"
df <- read.table(file, header = TRUE)
head(df)
EnhancedVolcano(df,
                lab = rownames(df),
                x = 'AR',
                y = 'qvalue',
                pCutoff = 0.05,
                labSize = 0.0,
                col=c('black', 'black', 'red3'),
                xlim = c(0, 1))
##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig5c.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
print(head(df))

# df$Region <- factor(df$Region,levels = df$Region,ordered = TRUE)
p <- ggplot(data=df, aes(x=Enrichment, y=Enrichment_p)) +
  # geom_line(color="#0072B2") +
  geom_point() + 
  geom_text(label=df$Trait,check_overlap=TRUE,hjust=1.2,nudge_x=0,family = "Arial",aes(colour=factor(color))) + 
  scale_colour_manual(values=c("#ff6666","#6699ff")) + 
  theme_bw() + labs(x="", y="-log10P (enrichment)") +
  # annotate("segment",x=0,xend=45,y=2,yend=2) + 
  theme(
    axis.title.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 14, family = "Arial", color = "black"),
    axis.text.x = element_text(family = "Arial", color = "black", size=12),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig5c.pdf", 
       p, dpi=600, width = 5.36, height = 5)
#########################################################################
library(forestplot)
# in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig5d.txt"
# df <- read.table(in_file, sep = "\t", header = TRUE)
# print(head(df))
#######
h2 <- 
  structure(list(
    mean = c(NA, 0.046560505, 0.054090487, 0.032869393, 0.040835625),
    lower = c(NA, 0.027110226, 0.030522386, 0.017005301, 0.011957489),
    upper = c(NA, 0.066010783, 0.077658588, 0.048733486,0.069713762)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -5L), 
    class = "data.frame")

tabletext <- cbind(
  c("Traits", "Schizophrenia", "BMI", "Years of education", "Dejerine-Sottas Disease"),
  c("Heritability(95% CI)", "0.05(0.03-0.07)", "0.05(0.03-0.08)", "0.03(0.02-0.05)", "0.04(0.01-0.07)")
)
forestplot(tabletext, 
           h2,new_page = TRUE,
           #is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
           # clip=c(0.1,2.5),
           # xlog=TRUE,
           #col=fpColors(box="royalblue",line="darkblue", summary="royalblue")
           )
##########################################################################
