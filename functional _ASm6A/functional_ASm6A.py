#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @author Shuo Cao, Yuhe Li

#This script aims to find variants which have same allelic direction in at least 90% of samples(located in m6A peak)
#Also 90% p-valus larger than 0.1 when performing fisher exact in sample pairs.
from collections import Counter
import os
import glob
import pandas as pd
from scipy.stats import fisher_exact
import statsmodels.stats.multitest
import numpy as np
#### sig + unsig (both)

ASm6A_dir = ""
Prop = ""
result_dir = ""

def select_functional_ASm6A(prop):
'''
The SNPs overlaped with m6A peak, did not filter by p-value, it should have same allelic direction in prop percent of samples.
In this article we use prop = 0.9
'''
    m6a_dict = {}
    bed_list = glob.glob("%s/*.bed" % ASm6A_dir)
    for bed in bed_list:
        if os.path.getsize(bed):
            sample = os.path.basename(bed)
            df = pd.read_table(bed, header=None)
            ## chr1:846078:ref>alt
            df['term'] = df.iloc[:,0] + ":" + df.iloc[:,2].astype(str) + ":" + df.iloc[:,3].str.split(";").str[0]
            df['mark'] = df.iloc[:,3].str.split(";").str[-1]
            term_list, mark_list = df['term'].tolist(), df['mark'].tolist()
            for i in range(len(term_list)):
                m6a_dict[term_list[i]] = m6a_dict.get(term_list[i], []) + [mark_list[i]]

    final_list = []
    for term, mark_list in m6a_dict.items():
        ## appear at least twice
        test_list = [x for x in mark_list if x != "unsig"]
        if len(test_list) >= 2: 
            count, new_line = Counter(mark_list), ""
            if (count['ref']/(count['ref']+count['alt']+count['unsig'])) >= prop: 
                new_line = [term.split(":")[0], int(term.split(":")[1])-1, 
                            int(term.split(":")[1]), "ref", term.split(":")[2], 
                           count['ref'], count['alt'], count['unsig']]
                final_list.append(new_line)
            elif (count['alt']/(count['ref']+count['alt']+count['unsig'])) >= prop: 
                new_line = [term.split(":")[0], int(term.split(":")[1])-1, 
                            int(term.split(":")[1]), "alt", term.split(":")[2], 
                           count['ref'], count['alt'], count['unsig']]
                final_list.append(new_line)
    df = pd.DataFrame(final_list)
    print(prop, len(df))
    return(df.drop_duplicates())

cm6a_df = select_functional_ASm6A(Prop)
functional_list = list(cm6a_df.iloc[:,0] + ":" + cm6a_df.iloc[:,2].astype(str) + ":" + cm6a_df.iloc[:,4])

def get_ratio(ASm6A_dir, functional_list) :
'''
get m6A level of ref allele and alt allele
refRPKM_ratio = ref_RPKM_IP/ref_RPKM_Input
altRPKM_ratio = alt_RPKM_IP/alt_RPKM_Input
'''
    file_dir = '%s/contained_m6A/merged_both/' % ASm6A_dir
    file_list = [os.path.basename(f) for f in glob.glob("%s/*.bed" % file_dir)]
    os.chdir(file_dir)
    ratio_dict = {}
    for file in file_list :
        file_df = pd.read_table(file, header=None)
        file_df.columns = ['chr','start','end','term','refRPKM_ratio','altRPKM_ratio','mark']
        term_list = list(file_df.iloc[:,0] + ":" + file_df.iloc[:,2].astype(str) + ":" + file_df.iloc[:,3].str.split(';').str[0])
        for c in functional_list :
            if c in term_list :
                idx = term_list.index(c)
                ratio_dict[c] = ratio_dict.get(c, []) + [file_df.loc[idx,['refRPKM_ratio','altRPKM_ratio']].tolist()]
    return ratio_dict
    
def fisher_each_samples(ratio_dict, result_dir, prop) :
'''
For all samples of one site, run fisher_exact and 90% of p-values should larger than 0.1
'''
    res_l = []
    for k,ratio_l in ratio_dict.items() :
        p_l = []
        
        for i in range(len(ratio_l)) :
            pl = []
            for j in range(i+1,len(ratio_l)) :
                p = fisher_exact([ratio_l[i], ratio_l[j]])[1]
                pl.append(p)
            zero = [0 for k in range(i+1)]
            p_l.append(zero + pl)
        p_arr = np.array(p_l)
        p_arr += p_arr.T - np.diag(p_arr.diagonal())
        
        counts = 0
        for i in range(len(ratio_l)) :
            idx_l = list(range(len(ratio_l)))
            idx_l.remove(i)
            sample_p = p_arr[i,idx_l]
            if sample_p[sample_p>0.1].shape[0]/(len(ratio_l)-1) >= 0.9 :
                counts += 1
            else :
                continue
        if counts == len(ratio_l) :
            res_l.append(k)
    with open('%s/functional_ASm6A_%s_fisher.bed' % (result_dir, prop), 'w') as out :
        for r in res_l :
            f = r.split(':')
            outline = '\t'.join([f[0], str(int(f[1])-1), f[1]])
            out.write(outline + '\n')

ratio_dict = get_ratio(ASm6A_dir, functional_list)
fisher_each_samples(ratio_dict, result_dir, prop)
