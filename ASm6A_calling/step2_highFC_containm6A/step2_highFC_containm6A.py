#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @author Shuo Cao

#step 2 
### SNPs with high foldchange and in m6A peak (remove m6Am)
### continue to filter
### 1. sig (qvalue < 0.05)
###    1.1 within m6A peak
###           1.1.1 ip/input >= 2 
### 2. unsig (qvalue < 0.05)
###     2.1 within m6A peak

ASm6A_dir = ""


def select_sig_and_unsig(in_file, data_type):
    df = pd.read_table(in_file)
    #############################
    if data_type == "sig":
        df = df[(df['qvalue']<0.05)] 
    elif data_type == "unsig":
        df = df[(df['qvalue']>=0.05)] 
    df = df.round(2)
    df['mark'] = np.where(df['altRPKM_ratio'] > df['refRPKM_ratio'], "alt", "ref")
    df['start'] = df['position'] - 1
    df['term'] = df['refAllele'] + ">" + df['altAllele'] + ";" + df['allelicRatio'].astype(str) + ";" + df['mark']
    df = df[['contig', 'start', 'position', 'term', 'refRPKM_ratio', 'altRPKM_ratio']]
    return df
    
os.chdir(ASm6A_dir)
sig_result_dir, unsig_result_dir = "sig/", "unsig/"
os.system("mkdir -p %s" % sig_result_dir)
os.system("mkdir -p %s" % unsig_result_dir)
file_list = glob.glob("*.txt")
for x in file_list:
    df_sig, df_unsig = select_sig_and_unsig(x, "sig"), select_sig_and_unsig(x, "unsig")
    print(os.path.basename(x), len(df_sig), len(df_unsig))
    res_sig = os.path.join(sig_result_dir, x.replace(".txt",".bed"))
    res_unsig = os.path.join(unsig_result_dir, x.replace(".txt",".bed"))
    df_sig.to_csv(res_sig, sep="\t", header=False, index=False)
    df_unsig.to_csv(res_unsig, sep="\t", header=False, index=False)
    
    
def overlap_with_m6A(in_dir, m6a_dir):
    os.chdir(in_dir)
    result_dir = "contained_m6A/"
    os.system("mkdir -p %s" % result_dir)

    #m6a_dir = "/Charles/project/ASm6A/peak_calling/merged_peak_MSPC/rm_m6Am/"

    bed_list = glob.glob("%s/*.bed" % in_dir)
    for asm6a in bed_list:
        prefix = os.path.basename(asm6a).split(".")[0]
        m6a = os.path.join(m6a_dir, "%s_m6A.bed" % prefix)
        res = os.path.join(result_dir, os.path.basename(asm6a))
        os.system("bedtools intersect -a %s -b %s -wa | sort -u > %s" % (asm6a, m6a, res))

        
        
overlap_with_m6A("%s/sig/" % ASm6A_dir)
overlap_with_m6A("%s/unsig/" % ASm6A_dir)


def select_highFC(in_dir):
    os.chdir(in_dir)
    result_dir = "highFC/"
    os.system("mkdir -p %s" % result_dir)
    bed_list = glob.glob("%s/*.bed" % in_dir)
    for asm6a in bed_list:
        df = pd.read_table(asm6a, header=None)
        df.columns = ['contig', 'start', 'position', 'term', 'refRPKM_ratio', 'altRPKM_ratio']
        df['mark'] = df['term'].str.split(";").str[2]
        #####################################################
        df = df[((df['mark']=="ref")&(df['refRPKM_ratio']>=2)) | ((df['mark']=="alt")&(df['altRPKM_ratio']>=2))]
        ####################################################
        res = os.path.join(result_dir, os.path.basename(asm6a))
        df.to_csv(res, sep="\t", header=False, index=False)
        
select_highFC("%s/sig/contained_m6A/" % ASm6A_dir)
