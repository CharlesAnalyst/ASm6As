library(ggplot2)
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

ggsave("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/VEF/pie.pdf", p, dpi=600, width = 30, height = 15, units = "cm")
##########################################
library(ggplot2)

in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig2a.txt"
df <- read.table(in_file, sep = "\t")
print(df)
df$V1 = factor(df$V1, levels = df$V1)

p <- ggplot(data=df, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.8)+
  scale_fill_manual(values=c("#D55E00", "#D55E00", "#D55E00", "#009E73"))+
  theme_bw() + guides(fill=FALSE)+ labs(x="", y="Log2 OR") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(angle = 60, hjust = 1, family = "Times", color = "black",face="bold", size=10),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig2a.pdf", 
       p, dpi=600, width = 3.6, height = 4.13)
########################################################
library(ggplot2)

in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig2c.txt"
df <- read.table(in_file, sep = "\t")
colnames(df) <- c("Function category", "RBP", "OR")
df<- df[seq(dim(df)[1],1),]
print(df)
df$RBP = factor(df$RBP, levels = df$RBP)
### 
#    n       reader     regul    trans       complex    
col_list = c("#F0E442","#0072B2","#009E73","#D55E00","#56B4E9")
# col_list = c("#009E73","#56B4E9","#009E73","#009E73","#56B4E9","#0072B2","#56B4E9","#009E73","#0072B2","#56B4E9",
#              "#0072B2","#56B4E9","#D55E00","#56B4E9","#56B4E9","#009E73","#56B4E9","#009E73","#56B4E9","#F0E442",
#              "#56B4E9","#009E73","#009E73","#009E73")
p <- ggplot(data=df, aes(x=RBP, y=OR, fill=`Function category`)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.8)+
  scale_fill_manual(values=col_list)+
  theme_bw() + coord_flip() + labs(x="", y="Log2 OR") + # + guides(fill=FALSE)
  theme(
    axis.title.x = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.y = element_text(size = 12, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(family = "Times", color = "black",face="bold", size=10),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig2c.pdf", 
       p, dpi=600, width = 5.91, height = 7.58)

###########################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/e_fig2a.txt"
df <- read.table(in_file, sep = "\t", header = TRUE)
# df<- df[seq(dim(df)[1],1),]
print(df)
df$RBP = factor(df$RBP, levels = df$RBP)

p <- ggplot(data=df, aes(x=RBP, y=Number, fill=Enriched)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.8) +
  scale_fill_manual(values=c("#0072B2","#D55E00"))+
  theme_bw() + labs(x="", y="Log2 OR") + # + guides(fill=FALSE)
  theme(
    axis.title.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.y = element_text(size = 12, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", color = "black",face="bold", size=10),
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
    axis.title.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.y = element_text(size = 12, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", color = "black",face="bold", size=10),
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
    axis.title.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.y = element_text(size = 12, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", color = "black",face="bold", size=10),
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
    axis.title.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", color = "black",face="bold", size=12),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
    # legend.title = element_blank()
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig3a.pdf", 
       p, dpi=600, width = 5.2, height = 9.1)

##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4a.txt"
df <- read.table(in_file, sep = "\t") # , header = TRUE
print(head(df))

df$V1 <- factor(df$V1,levels = df$V1,ordered = TRUE)
p <- ggplot(data=df, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() + labs(x="", y="ASm6A number") + guides(fill=FALSE) +
  theme(
    axis.title.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", color = "black",face="bold", size=12),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
    # legend.title = element_blank()
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig4a.pdf", 
       p, dpi=600, width = 5, height = 6)

##################################################################################
# in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4b.txt"
# df <- read.table(in_file, sep = "\t", header = TRUE)
# print(df)
# # df$RBP = factor(df$RBP, levels = df$RBP)
# ### 
# #    n       reader     regul    trans       complex    
# col_list = c("#F0E442","#0072B2","#009E73","#D55E00","#56B4E9")
# # col_list = c("#009E73","#56B4E9","#009E73","#009E73","#56B4E9","#0072B2","#56B4E9","#009E73","#0072B2","#56B4E9",
# #              "#0072B2","#56B4E9","#D55E00","#56B4E9","#56B4E9","#009E73","#56B4E9","#009E73","#56B4E9","#F0E442",
# #              "#56B4E9","#009E73","#009E73","#009E73")
# p <- ggplot(data=df, aes(x=Description, y=X.log10p.adjust, fill=tissue)) +
#   geom_bar(stat="identity", position=position_dodge(), width = 0.8)+
#   scale_fill_manual(values=col_list)+
#   theme_bw() + coord_flip() + labs(x="", y="Log2 OR") + # + guides(fill=FALSE)
#   theme(
#     axis.title.x = element_text(size = 14, family = "Times", color = "black",face="bold"),
#     axis.text.y = element_text(size = 12, family = "Times", color = "black",face="bold"),
#     axis.text.x = element_text(family = "Times", color = "black",face="bold", size=10),
#     panel.border = element_blank(),
#     axis.line = element_line(colour = "black")
#   )
# p
# ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig4b.pdf", 
#        p, dpi=600, width = 5.91, height = 7.58)
##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4c.txt"
df <- read.table(in_file, sep = "\t") # , header = TRUE
print(head(df))

df$V1 <- factor(df$V1,levels = df$V1,ordered = TRUE)
p <- ggplot(data=df, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() + labs(x="", y="ASm6A number") + guides(fill=FALSE) +
  theme(
    axis.title.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", color = "black",face="bold", size=12),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
    # legend.title = element_blank()
  )
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/fig4c.pdf", 
       p, dpi=600, width = 5.2, height = 6.2)

##################################################################################
in_file = "/home/galaxy/project/alleleSpecific_analysis/plots/fig4d.txt"
df <- read.table(in_file, sep = "\t") # , header = TRUE
print(head(df))

# df$V1 <- factor(df$V1,levels = df$V1,ordered = TRUE)
p <- ggplot(data=df, aes(x=V2, y=V3, fill=V1)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() + labs(x="", y="ASm6A number") + guides(fill=FALSE) +
  theme(
    axis.title.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", color = "black",face="bold", size=12),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
    # legend.title = element_blank()
  ) + scale_x_discrete(limits=c("promoter-TSS", "5UTR", "exon", "intron", "TTS", "3UTR", "Intergenic", "non-coding"))
p
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
    axis.title.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.y = element_text(size = 14, family = "Times", color = "black",face="bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", color = "black",face="bold", size=12),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
    # legend.title = element_blank()
  ) + scale_x_discrete(limits=c("protein_coding", "antisense", "pseudogene", "processed_transcript", "sense_overlapping", "lincRNA", "3prime_overlapping_ncrna", "sense_intronic", "TR_C_gene"))
p
ggsave("/home/galaxy/project/alleleSpecific_analysis/plots/e_fig4a.pdf", 
       p, dpi=600, width = 8.3, height = 4.2)

###########################################################################
library(sunburstR)
library(htmltools)
library(d3r)
dat <- data.frame(
  level1 = rep(c("a", "b"), each=3),
  level2 = paste0(rep(c("a", "b"), each=3), 1:3),
  size = c(10,5,2,3,8,6),
  stringsAsFactors = FALSE
)

knitr::kable(dat)
tree <- d3_nest(dat, value_cols = "size")
tree
sb1 <- sunburst(tree, width="100%", height=400)
sb2 <- sunburst(
  tree,
  legend = FALSE,
  width = "100%",
  height = 400
)

# do side-by-side for comparison
div(
  style="display: flex; align-items:center;",
  div(style="width:50%; border:1px solid #ccc;", sb1),
  div(style="width:50%; border:1px solid #ccc;", sb2)
)
sb3 <- sund2b(tree, width="100%")


div(
  style="display: flex; align-items:center;",
  sb3
)

