import os
import re
import glob
import numpy as np
import pandas as pd
from scipy import stats
from collections import defaultdict, Counter


snp_list = glob.glob(os.path.join(hap1_hap2_dir,'*'))

###  fisher test 
with open('hap1_hap2_diff.fisher.txt','w') as f:
    for sub_dir in snp_list:
        hap1_file = os.path.join(sub_dir,'hap1.totalm6A.FDR.out')
        hap2_file = os.path.join(sub_dir,'hap2.totalm6A.FDR.out')
        hap1_table = pd.read_csv(hap1_file,sep='\t',header=None)
        hap1_table.columns = ['Chr','Sites','Strand','Gene','CR','AGcov','Acov','Genecov','Ratio','Pvalue','P_adjust']
        hap2_table = pd.read_csv(hap2_file,sep='\t',header=None)
        hap2_table.columns = ['Chr','Sites','Strand','Gene','CR','AGcov','Acov','Genecov','Ratio','Pvalue','P_adjust']


        hap1_A_dict = dict(zip(hap1_table['Sites'],hap1_table['Acov']))
        hap2_A_dict = dict(zip(hap2_table['Sites'],hap2_table['Acov']))
        hap1_AG_dict = dict(zip(hap1_table['Sites'],hap1_table['AGcov']))
        hap2_AG_dict = dict(zip(hap2_table['Sites'],hap2_table['AGcov']))

        common_m6a = set(hap1_A_dict.keys()).intersection(set(hap2_A_dict.keys()))

        for m6a in common_m6a:
            statistic,pvalue = stats.fisher_exact([[hap1_A_dict[m6a],hap2_A_dict[m6a]],[hap1_AG_dict[m6a],hap2_AG_dict[m6a]]])
            if pvalue < 0.05 and statistic > 1:
                f.write('{}\t{}\t{}\n'.format(os.path.basename(sub_dir),m6a,'ref'))
            elif pvalue < 0.05 and statistic < 1:
                f.write('{}\t{}\t{}\n'.format(os.path.basename(sub_dir),m6a,'alt'))
            else:
                f.write('{}\t{}\t{}\n'.format(os.path.basename(sub_dir),m6a,'unsig'))


###  ratio diff > 0.1
with open('hap1_hap2_diff.0.1.txt','w') as f:
    for sub_dir in snp_list:
        hap1_file = os.path.join(sub_dir,'hap1.totalm6A.FDR.out')
        hap2_file = os.path.join(sub_dir,'hap2.totalm6A.FDR.out')
        hap1_table = pd.read_csv(hap1_file,sep='\t',header=None)
        hap1_table.columns = ['Chr','Sites','Strand','Gene','CR','AGcov','Acov','Genecov','Ratio','Pvalue','P_adjust']
        hap2_table = pd.read_csv(hap2_file,sep='\t',header=None)
        hap2_table.columns = ['Chr','Sites','Strand','Gene','CR','AGcov','Acov','Genecov','Ratio','Pvalue','P_adjust']

        hap1_dict = dict(zip(hap1_table['Sites'],hap1_table['Ratio']))
        hap2_dict = dict(zip(hap2_table['Sites'],hap2_table['Ratio']))

        common_m6a = set(hap1_dict.keys()).intersection(set(hap2_dict.keys()))
        for m6a in common_m6a:
            if hap1_dict[m6a] - hap2_dict[m6a] > 0.1:
                f.write('{}\t{}\t{}\t{}\n'.format(os.path.basename(sub_dir),m6a,hap1_dict[m6a] - hap2_dict[m6a],'ref'))
            elif hap2_dict[m6a] - hap1_dict[m6a] > 0.1:
                f.write('{}\t{}\t{}\t{}\n'.format(os.path.basename(sub_dir),m6a,hap1_dict[m6a] - hap2_dict[m6a],'alt'))
            else:
                f.write('{}\t{}\t{}\t{}\n'.format(os.path.basename(sub_dir),m6a,hap1_dict[m6a] - hap2_dict[m6a],'unsig'))
