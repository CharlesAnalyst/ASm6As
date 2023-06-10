# rm(list=ls())
library('ggplot2')
setwd('/home/galaxy/project/m6AQTL/data/data_for_plot')
data = read.table('tissue_m6a_count.txt',sep = '\t' , quote = NULL,header = F)
data<- data[order(data$V2,decreasing = T),]
data$color <- rep(c('1','2','3','4','5','6'),times = c( count(data[data$V2>5000,]),count(data[data$V2 <= 5000 & data$V2 >4000,]),count(data[data$V2 <= 4000 & data$V2 >3000,]),count(data[data$V2 <= 3000 & data$V2 >2000,]),count(data[data$V2 <= 2000 & data$V2 >1000,]),count(data[data$V2<=1000,])))
data<- data[order(data$V2),]
data$V1 <- factor(data$V1,levels = data$V1)
ggplot(data,aes(x = V1 ,y = V2,fill = color)) +geom_bar(stat = 'identity') +scale_fill_brewer(palette = 'Spectral',direction = -1)+
  ggtitle('m6A number')+geom_text(aes(label=V2),family = 'Times',hjust=0)+ 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle =45,hjust = 1,vjust = 1),axis.text=element_text(family = 'Times',size = 20,color = 'black',face = 'bold'),legend.title=element_blank(),axis.title = element_blank())+
  scale_y_continuous(limits = c(0,10000),expand=c(0.01,0),breaks = c(2000,4000,6000,8000))+guides(fill=FALSE)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  coord_flip()

rm(list=ls())
library('ggplot2')
setwd('/home/lzg/1906/m6AQTL/script')
gsmdata = read.table('tissue_gsm_count.txt',row.names = 1,sep = '\t' , quote = NULL,header = F)
data = read.table('tissue_m6a_count.txt',row.names = 1,sep = '\t' , quote = NULL,header = F)
all_data<-cbind(data,gsmdata[rownames(data),])
all_data<- all_data[order(all_data$V2,decreasing = T),]
all_data$color <- rep(c('1','2','3','4','5','6'),times = c( count(all_data[all_data$V2>5000,]),count(all_data[all_data$V2 <= 5000 & all_data$V2 >4000,]),count(all_data[all_data$V2 <= 4000 & all_data$V2 >3000,]),count(all_data[all_data$V2 <= 3000 & all_data$V2 >2000,]),count(all_data[all_data$V2 <= 2000 & all_data$V2 >1000,]),count(all_data[all_data$V2<=1000,])))
all_data<- all_data[order(all_data$V2),]
all_data$V1 <- factor(rownames(all_data),levels = rownames(all_data))
colnames(all_data)<-c('m6a_count','gsm_count','color','tissue')

ggplot(all_data,aes(x = tissue ,y = gsm_count,fill = color)) +geom_bar(stat = 'identity') +scale_fill_brewer(palette = 'Spectral')+
  ggtitle('GSM number')+ geom_text(aes(label=gsm_count),size = 9,family = 'Times',hjust=-0.5)+
  theme(plot.title = element_text(hjust = 0.3,family = 'Times',size = 30,face = 'bold'),axis.text=element_text(family = 'Times',size = 25,color = 'black',face = 'bold'),legend.title=element_blank(),axis.title = element_blank())+
  scale_y_continuous(limits = c(0,20),expand=c(0.01,0),breaks = c(0,5,10,15))+guides(fill=FALSE)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  coord_flip()

ggplot(all_data,aes(x = tissue ,y = m6a_count,fill = color)) +geom_bar(stat = 'identity') +scale_fill_brewer(palette = 'Spectral')+
  ggtitle('m6A number')+ geom_text(aes(label=m6a_count),size = 9,family = 'Times',hjust=-0.1)+
  theme(plot.title = element_text(hjust = 0.3,family = 'Times',size = 30,face = 'bold'),axis.text=element_text(family = 'Times',size = 25,color = 'black',face = 'bold'),legend.title=element_blank(),axis.title = element_blank())+
  scale_y_continuous(limits = c(0,11000),expand=c(0.01,0),breaks = c(0,2000,4000,6000,8000))+guides(fill=FALSE)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  coord_flip()

