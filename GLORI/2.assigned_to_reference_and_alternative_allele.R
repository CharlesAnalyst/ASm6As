library(asSeq)

work_dir = 'work_dir'
sub_dir_list = list.dirs(work_dir)
#bam = 'merged.sorted.bam'

for (each_dir in sub_dir_list){
  snplist = paste0(each_dir,'/snplist')
  bam = paste0(each_dir,'/snp.bam')
  extractAsReads(bam,snplist,outputTag = paste0(each_dir,'/var'))