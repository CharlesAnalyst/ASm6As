setwd("/home/galaxy/project/pausing_index/macs_result/for_plot_wt")
data=read.table('data_for_plot.wt',header=F,sep="\t")
peak_freq=c(as.numeric(as.vector(data[which(data[,8]=="Genes_TSS3k"),5])),as.numeric(as.vector(data[which(data[,8]=="Genes"),5]))+100,as.numeric(as.vector(data[which(data[,8]=="Genes_TTS3k"),5]))+200)
plot(density(peak_freq),xlim=c(0,300),lwd=3,ylab="Peak density",xlab=NA,main=NA,col="red",xaxt="n",type="l")
abline(v=100,lty=2,col="black",lwd=1)
abline(v=200,lty=2,col="black",lwd=1)
text(x=c(0,100,200,300),y=-0.001,pos=2,labels=c("-3k","TSS",'TTS',"+3k"),xpd=TRUE,font=2)

