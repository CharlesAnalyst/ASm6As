library(ciRcus)
setwd("/home/galaxy/project/circRNA_analysis/find_circ/total_files")
gtf2sqlite(assembly = "mm10", db.file  = "mouse_mm10_txdb.sqlite")
annot.list <- loadAnnotation("mouse_mm10_txdb.sqlite")
cdata <- data.frame(sample = c("ko", "ko2", "wt", "wt2"), filename=list.files("/home/galaxy/project/circRNA_analysis/find_circ/total_files/find_circ_bed"), pattern=".bed", full.name=TRUE)
cdata
circs.se <- summarizeCircs(colData=cdata, wobble=1, keepCols=1:12)
circs.se <- annotateCircs(circs.se, annot.list=annot.list)
circs.dt <- resTable(circs.se)
write.table(circs.dt, "CIRCUS_result.txt", sep = "\t", quote = FALSE, row.names = FALSE)
## plot
histogram(circs.se, 0.5)
annotPie(circs.se, 0.02)
uniqReadsQC(circs.se, "all")
