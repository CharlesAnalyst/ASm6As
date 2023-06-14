#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @author Shuo Cao

### filter order
### remove SNPs which don't meet minimum needs
### Hypothesis test

from scipy import stats
from multiprocessing import Pool
from statsmodels.stats.multitest import fdrcorrection

count_dir = ""
result_dir = ""
map_file = ""
summary_file = ""

os.chdir(count_dir)
os.system("mkdir -p %s" % result_dir)

#step 1 limit IP total number and Hypothesis test
### 1. only limit IP total number 
###    1.1 input (refCount, altCount >=2 and refCount+altCount >=10)
###    1.2 ip (refCount+altCount >=10)

def generate_totalNum_dict(map_file, summary_file):
'''
generate dict which keys are prefix of IP files and which values are lists contain number of total reads in IP files, prefix of Input files and number of total reads in Input files
'''
    df_map = pd.read_table(map_file)
    map_dict = dict(zip(df_map['IP'], df_map['Input']))
    df_stat = pd.read_table(summary_file, header=None)
    number_dict = dict(zip(df_stat.iloc[:,0].tolist(), df_stat.iloc[:,1].tolist()))
    totalNum_dict = {}
    for sample in number_dict.keys():
        if sample in map_dict:
            totalNum_dict[sample] = [number_dict[sample], map_dict[sample], number_dict[map_dict[sample]]]
    return totalNum_dict

def preprocess_df(ip_prefix):
'''
RPKM normalization
'''
    input_prefix = totalNum_dict[ip_prefix][1]
    ip_file, input_file = "ip_count/%s.readcounts.txt" % ip_prefix, "input_count/%s.readcounts.txt" % input_prefix
    df_ip, df_input = pd.read_table(ip_file), pd.read_table(input_file)
    df_results = []
    totalNum_list = [totalNum_dict[ip_prefix][0], totalNum_dict[ip_prefix][2]]
    treat_list, df_list = ['ip', 'input'], [df_ip, df_input]
    for i in range(len(df_list)):
        df, treat, total_num = df_list[i], treat_list[i], totalNum_list[i]
        df = df[['contig', 'position', 'refAllele', 'altAllele', 'refCount', 'altCount']]
        #######################################################
        if i == 0: # ip sample: only restrict total number
            df = df[(df['refCount'] + df['altCount'])>=10]
        else:
            df = df[(df['refCount']>=2) & (df['altCount']>=2)]
            df = df[(df['refCount'] + df['altCount'])>=10]
        #####################################################################
        df['refRPKM_%s' % treat] = (df['refCount'] / total_num) * 1000000000
        df['altRPKM_%s' % treat] = (df['altCount'] / total_num) * 1000000000
        df_results.append(df)
        
    return df_results[0], df_results[1]

def fisher_exact_test_each(ip_prefix):
'''
Hypothesis test and FDR correction
'''
    result_file = os.path.join(result_dir, ip_prefix+".txt")
    df_ip, df_input = preprocess_df(ip_prefix)
    df = df_ip.merge(df_input, on=['contig', 'position', 'refAllele', 'altAllele'], how='left').dropna(how="any")
    print(ip_prefix, len(df))
    df['refRPKM_ratio'] = df['refRPKM_ip'] / df['refRPKM_input']
    df['altRPKM_ratio'] = df['altRPKM_ip'] / df['altRPKM_input']
    df['allelicRatio'] = df['refRPKM_ratio'] / (df['altRPKM_ratio'] + df['refRPKM_ratio'])
    pvalue_list , odds_list = [], []
    for i,j in df.iterrows():
        a,b = j['refRPKM_ip']+0.5, j['altRPKM_ip']+0.5
        c,d = j['refRPKM_input']+0.5, j['altRPKM_input']+0.5
        oddsratio,pvalue=stats.fisher_exact([[a,b], [c,d]])
        pvalue_list.append(pvalue)
        odds_list.append(oddsratio)
    qvalue_list = list(fdrcorrection(pvalue_list, method="indep")[1]) # fdr_bh
    df['pvalue'],df['oddsratio'],df['qvalue'] = pvalue_list,odds_list,qvalue_list
    df_sub=df[['contig','position','refAllele','altAllele','refRPKM_ratio','altRPKM_ratio', 'allelicRatio', 'pvalue', 'oddsratio','qvalue']]
    df_sub.to_csv(result_file, sep="\t", index=False)

totalNum_dict = generate_totalNum_dict(map_file, summary_file)
prefix_list = totalNum_dict.keys()

pool = Pool(processes=30)
for ip_prefix in prefix_list:
    pool.apply_async(fisher_exact_test_each, (ip_prefix, ))
pool.close()
pool.join()


