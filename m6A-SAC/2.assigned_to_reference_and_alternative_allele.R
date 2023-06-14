library(asSeq)

work_dir = 'work_dir'
sub_dir_list = list.dirs(work_dir)

for (each_dir in sub_dir_list){
  snplist = paste0(each_dir,'/snplist')
  original_bam_dir = paste0(each_dir,'/original_bam_dir')
  
  input_plus_snp_bam = paste0(original_bam_dir,'/input_plus_snp.bam')
  input_minus_snp_bam = paste0(original_bam_dir,'/input_minus_snp.bam')
  ftom_rep1_plus_snp_bam = paste0(original_bam_dir,'/ftom_rep1_plus_snp.bam')
  ftom_rep1_minus_snp_bam = paste0(original_bam_dir,'/ftom_rep1_minus_snp.bam')
  ftom_rep2_plus_snp_bam = paste0(original_bam_dir,'/ftom_rep2_plus_snp.bam')
  ftom_rep2_minus_snp_bam = paste0(original_bam_dir,'/ftom_rep2_minus_snp.bam')
  
  out_bam_dir = paste0(each_dir,'/allele_bam_dir')
  extractAsReads(input_plus_snp_bam,snplist,outputTag = paste0(out_bam_dir,'/input_plus_var'))
  extractAsReads(input_minus_snp_bam,snplist,outputTag = paste0(out_bam_dir,'/input_minus_var'))
  extractAsReads(ftom_rep1_plus_snp_bam,snplist,outputTag = paste0(out_bam_dir,'/ftom_rep1_plus_var'))
  extractAsReads(ftom_rep1_minus_snp_bam,snplist,outputTag = paste0(out_bam_dir,'/ftom_rep1_minus_var'))
  extractAsReads(ftom_rep2_plus_snp_bam,snplist,outputTag = paste0(out_bam_dir,'/ftom_rep2_plus_var'))
  extractAsReads(ftom_rep2_minus_snp_bam,snplist,outputTag = paste0(out_bam_dir,'/ftom_rep2_minus_var'))
}