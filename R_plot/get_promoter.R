library("TxDb.Hsapiens.UCSC.hg38.knownGene")

setwd("/data5/galaxy/project/data/promoter")
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- BSgenome.Hsapiens.NCBI.GRCh38
genes_txdb <- ensgene(txdb)
data(TSS.human.GRCh38)
promoters_txdb <- promoters(genes_txdb)
#By default, the promoters function will fetch the 2000 nucleotides before the transcription start site (TSS) and the 200 nucleotides after the TSS.
promoters_txdb
#DataFrame()