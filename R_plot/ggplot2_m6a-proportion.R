library(ggplot2)
setwd("/data5/galaxy/project/CpG_proporation-m6aGene")
data = read.table("plot_data.txt", sep = "\t", header = TRUE)
head(data)
p <- ggplot(data, aes(x = CpG_type, y = m6a.propor, fill = CpG_type)) + geom_boxplot() + scale_x_discrete(limits=c("low CpG", "high CpG"))
p1 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line())
p2 <- p1 + scale_fill_manual(values = c("#87CEEB", "#EE7600"), limits=c("low CpG", "high CpG"))
p3 <- p2 + ylab("m6a proportion") + xlab("") + theme(axis.title.y = element_text(family="Times",face="bold", size = 18))
p4 <- p3 + theme(axis.text.x=element_text(family="Times", face="bold", color="black", size=rel(1.4))) + theme(axis.text.y=element_text(family="Times", face="bold", color="black", size=rel(1.2)))
p5 <- p4 + ylim(0.0, max(data$m6a.propor)) + labs(fill="")
ggsave("ggplot2_boxplot.pdf", width = 4, height = 4)

