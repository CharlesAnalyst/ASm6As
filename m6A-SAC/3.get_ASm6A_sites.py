import os
import re
import glob
import numpy as np
import pandas as pd
from scipy import stats
from collections import defaultdict, Counter


###  fisher test 
snp_list = os.listdir(hap1_hap2_dir)
each_snp_dir = 'each_snp_dir'
with open('ftom_rep1_input.fisher.txt','w') as f:
    for snp in snp_list:
        # print(snp)
        snp_chrom,snp_start,snp_end = snp.split('_')
        work_dir = os.path.join(each_snp_dir,snp)
        identification_dir = os.path.join(work_dir,'identification_dir')
        ftom_rep1_input_hap1_site = os.path.join(identification_dir,'ftom_rep1_input_hap1_DRACH.site')
        ftom_rep1_input_hap2_site = os.path.join(identification_dir,'ftom_rep1_input_hap2_DRACH.site')
        if os.path.exists(ftom_rep1_input_hap1_site) and os.path.exists(ftom_rep1_input_hap2_site):
            if os.path.getsize(ftom_rep1_input_hap1_site) > 0 and os.path.getsize(ftom_rep1_input_hap2_site) > 0:
                hap1_table = pd.read_csv(ftom_rep1_input_hap1_site,sep='\t',header=None)
                hap1_table.columns = ['Chr','Start','End','base','pos','strand','level','motif','input_ref_reads','input_alt_reads','input_var_freq','ftom_ref_reads','ftom_alt_reads','ftom_var_freq','p1','pvalue']
                hap2_table = pd.read_csv(ftom_rep1_input_hap2_site,sep='\t',header=None)
                hap2_table.columns = ['Chr','Start','End','base','pos','strand','level','motif','input_ref_reads','input_alt_reads','input_var_freq','ftom_ref_reads','ftom_alt_reads','ftom_var_freq','p1','pvalue']
                
                common_sites = list(set(hap1_table['pos'].tolist()).intersection(set(hap2_table['pos'].tolist())))
                for common_site in common_sites:
                    hap1_level = hap1_table.loc[hap1_table['pos'] == common_site,'level'].values[0]
                    hap2_level = hap2_table.loc[hap2_table['pos'] == common_site,'level'].values[0]
                    hap1_ftom_ref_reads = hap1_table.loc[hap1_table['pos'] == common_site,'ftom_ref_reads'].values[0]
                    hap1_ftom_alt_reads = hap1_table.loc[hap1_table['pos'] == common_site,'ftom_alt_reads'].values[0]
                    hap2_ftom_ref_reads = hap2_table.loc[hap2_table['pos'] == common_site,'ftom_ref_reads'].values[0]
                    hap2_ftom_alt_reads = hap2_table.loc[hap2_table['pos'] == common_site,'ftom_alt_reads'].values[0]

                    statistic,pvalue = stats.fisher_exact([[hap1_ftom_alt_reads,hap2_ftom_alt_reads],[(hap1_ftom_ref_reads+hap1_ftom_alt_reads),(hap2_ftom_ref_reads+hap2_ftom_alt_reads)]])

                    if pvalue < 0.05 and statistic > 1:
                        chrom = hap1_table.loc[hap1_table['pos'] == common_site,'Chr'].values[0]
                        start = hap1_table.loc[hap1_table['pos'] == common_site,'Start'].values[0]
                        end = hap1_table.loc[hap1_table['pos'] == common_site,'End'].values[0]
                        ftom_var_freq = hap1_table.loc[hap1_table['pos'] == common_site,'ftom_var_freq'].values[0]
                        input_var_freq = hap1_table.loc[hap1_table['pos'] == common_site,'input_var_freq'].values[0]
                        pvalue = hap1_table.loc[hap1_table['pos'] == common_site,'pvalue'].values[0]
                        ftom_alt_reads = hap1_table.loc[hap1_table['pos'] == common_site,'ftom_alt_reads'].values[0]
                        if (ftom_var_freq - input_var_freq > 5 and pvalue < 0.05) or (ftom_var_freq - input_var_freq > 10 and ftom_alt_reads >= 5):
                            f.write(snp_chrom+'\t'+snp_start+'\t'+snp_end+'\t'+common_site+'\t'+str(hap1_level - hap2_level)+'\t'+'ref'+'\n')
                    elif pvalue < 0.05 and statistic < 1:
                        chrom = hap1_table.loc[hap1_table['pos'] == common_site,'Chr'].values[0]
                        start = hap1_table.loc[hap1_table['pos'] == common_site,'Start'].values[0]
                        end = hap1_table.loc[hap1_table['pos'] == common_site,'End'].values[0]
                        ftom_var_freq = hap2_table.loc[hap2_table['pos'] == common_site,'ftom_var_freq'].values[0]
                        input_var_freq = hap2_table.loc[hap2_table['pos'] == common_site,'input_var_freq'].values[0]
                        pvalue = hap2_table.loc[hap2_table['pos'] == common_site,'pvalue'].values[0]
                        ftom_alt_reads = hap2_table.loc[hap2_table['pos'] == common_site,'ftom_alt_reads'].values[0]
                        if (ftom_var_freq - input_var_freq > 5 and pvalue < 0.05) or (ftom_var_freq - input_var_freq > 10 and ftom_alt_reads >= 5):
                            f.write(snp_chrom+'\t'+snp_start+'\t'+snp_end+'\t'+common_site+'\t'+str(hap1_level - hap2_level)+'\t'+'alt'+'\n')


###  ratio diff > 0.1
snp_list = os.listdir(hap1_hap2_dir)
each_snp_dir = 'each_snp_dir'
with open('ftom_rep1_input.0.1.txt','w') as f:
    for snp in snp_list:
        # print(snp)
        snp_chrom,snp_start,snp_end = snp.split('_')
        work_dir = os.path.join(each_snp_dir,snp)
        identification_dir = os.path.join(work_dir,'identification_dir')
        ftom_rep1_input_hap1_site = os.path.join(identification_dir,'ftom_rep1_input_hap1_DRACH.site')
        ftom_rep1_input_hap2_site = os.path.join(identification_dir,'ftom_rep1_input_hap2_DRACH.site')
        if os.path.exists(ftom_rep1_input_hap1_site) and os.path.exists(ftom_rep1_input_hap2_site):
            if os.path.getsize(ftom_rep1_input_hap1_site) > 0 and os.path.getsize(ftom_rep1_input_hap2_site) > 0:
                # print(snp)
                hap1_table = pd.read_csv(ftom_rep1_input_hap1_site,sep='\t',header=None)
                hap1_table.columns = ['Chr','Start','End','base','pos','strand','level','motif','input_ref_reads','input_alt_reads','input_var_freq','ftom_ref_reads','ftom_alt_reads','ftom_var_freq','p1','pvalue']
                hap2_table = pd.read_csv(ftom_rep1_input_hap2_site,sep='\t',header=None)
                hap2_table.columns = ['Chr','Start','End','base','pos','strand','level','motif','input_ref_reads','input_alt_reads','input_var_freq','ftom_ref_reads','ftom_alt_reads','ftom_var_freq','p1','pvalue']
                
                common_sites = list(set(hap1_table['pos'].tolist()).intersection(set(hap2_table['pos'].tolist())))
                for common_site in common_sites:
                    hap1_level = hap1_table.loc[hap1_table['pos'] == common_site,'level'].values[0]
                    hap2_level = hap2_table.loc[hap2_table['pos'] == common_site,'level'].values[0]
                    if hap1_level - hap2_level > 10:
                        chrom = hap1_table.loc[hap1_table['pos'] == common_site,'Chr'].values[0]
                        start = hap1_table.loc[hap1_table['pos'] == common_site,'Start'].values[0]
                        end = hap1_table.loc[hap1_table['pos'] == common_site,'End'].values[0]
                        ftom_var_freq = hap1_table.loc[hap1_table['pos'] == common_site,'ftom_var_freq'].values[0]
                        input_var_freq = hap1_table.loc[hap1_table['pos'] == common_site,'input_var_freq'].values[0]
                        pvalue = hap1_table.loc[hap1_table['pos'] == common_site,'pvalue'].values[0]
                        ftom_alt_reads = hap1_table.loc[hap1_table['pos'] == common_site,'ftom_alt_reads'].values[0]
                        if (ftom_var_freq - input_var_freq > 5 and pvalue < 0.05) or (ftom_var_freq - input_var_freq > 10 and ftom_alt_reads >= 5):
                            f.write(snp_chrom+'\t'+snp_start+'\t'+snp_end+'\t'+common_site+'\t'+str(hap1_level - hap2_level)+'\t'+'ref'+'\n')
                    elif hap2_level - hap1_level > 10:
                        chrom = hap1_table.loc[hap1_table['pos'] == common_site,'Chr'].values[0]
                        start = hap1_table.loc[hap1_table['pos'] == common_site,'Start'].values[0]
                        end = hap1_table.loc[hap1_table['pos'] == common_site,'End'].values[0]
                        ftom_var_freq = hap2_table.loc[hap2_table['pos'] == common_site,'ftom_var_freq'].values[0]
                        input_var_freq = hap2_table.loc[hap2_table['pos'] == common_site,'input_var_freq'].values[0]
                        pvalue = hap2_table.loc[hap2_table['pos'] == common_site,'pvalue'].values[0]
                        ftom_alt_reads = hap2_table.loc[hap2_table['pos'] == common_site,'ftom_alt_reads'].values[0]
                        if (ftom_var_freq - input_var_freq > 5 and pvalue < 0.05) or (ftom_var_freq - input_var_freq > 10 and ftom_alt_reads >= 5):
                            f.write(snp_chrom+'\t'+snp_start+'\t'+snp_end+'\t'+common_site+'\t'+str(hap1_level - hap2_level)+'\t'+'alt'+'\n')
