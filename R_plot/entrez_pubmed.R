# BiocManager::install("easyPubMed")
# library(easyPubMed)
# library(dplyr)
# library(kableExtra)


########################################################################################################
library(rentrez)
# term <- '(("sQTL") OR ("sQTLs") OR ("splicing quantitative trait loci")) AND 2008:2019[PDAT] AND human[ORGN]'
# term <- '("sQTL") OR ("sQTLs") AND human[ORGN]'
term <- '("pQTL") OR ("pQTLs") AND human[ORGN]'
# term <- 'Homo[ORGN] AND 2008:2019[PDAT] AND (("meQTL") OR ("meQTLs"))'
# query <- paste(term, "AND (", year, "[PDAT])")
# print(query)
res <- entrez_search(db = "pubmed", term = term, retmax=10000)
res$ids
res$count

year <- 2008:2019
papers <- sapply(year, search_year)
search_year <- function(year, term){
  query <- paste(term, "AND (", year, "[PDAT])")
  entrez_search(db="pubmed", term=query, retmax=0)$count
}
