import pandas as pd
from collections import OrderedDict
import numpy as np
from matplotlib_venn import venn2,venn2_circles
from matplotlib import pyplot as plt
import os
import argparse
plt.rcParams['pdf.fonttype'] = 42

def read_fasta_to_dic(filename):
    """
    function used to parser small fasta
    still effective for genome level file
    """
    fa_dic = OrderedDict()

    with open(filename, "r") as f:
        for n, line in enumerate(f.readlines()):
            if line.startswith(">"):
                if n > 0:
                    fa_dic[full_name] = "".join(seq_l)  # store previous one

                full_name = line.strip().replace(">", "")
                short_name = full_name.split(" ")[0]
                seq_l = []
            else:  # collect the seq lines
                if len(line) > 8:  # min for fasta file is usually larger than 8
                    seq_line1 = line.strip()
                    seq_l.append(seq_line1)

        fa_dic[short_name] = "".join(seq_l)  # store the last one
    return fa_dic
def extract_tombo(ab_path,fasta_path):
    tombo_result_minus_path = ab_path + ".statistic.minus.wig"
    tombo_result_plus_path = ab_path + ".statistic.plus.wig"
    tombo_difference_minus_path = ab_path + ".difference.minus.wig"
    tombo_difference_plus_path = ab_path + ".difference.plus.wig"
    ref_dict = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    def read_tombo_results(path):
        result_table=[]
        f = open(path,'r')
        chrome = None
        for line in f:
            if line[0:5]=='track':
                continue
            if 'chrom' in line:
                chrome = line.split(' ')[1].split('=')[1]
            else:
                static_list = line[:-1].split(' ')
                result_table.append([chrome,static_list[0],static_list[1]])
            if chrome is None:
                raise RuntimeError("Error")

        df = pd.DataFrame(result_table)
        # df = df[df[1]>=0.2]
        df[1] = df[1].astype(int)
        df[2]=df[2].astype(float)
        df.insert(loc=2, column="A1", value=df[1].values + 1)
        df.insert(loc=3, column="A2", value=np.repeat([["."]], df.shape[0]))
        df.insert(loc=4, column="A3", value=np.repeat([["."]], df.shape[0]))
        # sns.displot(df, x=1)
        # plt.show()
        if path.find("plus") != -1:
            # df = df[df[2] == "A"]
            df.insert(loc=5, column="A4", value=np.repeat([["+"]], df.shape[0]))
        else:
            # df = df[df[2] == "T"]
            df.insert(loc=5, column="A4", value=np.repeat([["-"]], df.shape[0]))
        df.columns=[0,1,2,3,4,5,6]
        return df

    plus_df = read_tombo_results(tombo_difference_plus_path)
    try:
        minus_df = read_tombo_results(tombo_difference_minus_path)
    except Exception:
        minus_df = None
    df_new = pd.concat([plus_df, minus_df], ignore_index=True)

    plus_df = read_tombo_results(tombo_result_plus_path)
    try:
        minus_df = read_tombo_results(tombo_result_minus_path)
    except Exception:
        minus_df = None
    df = pd.concat([plus_df, minus_df], ignore_index=True)
    all_df = pd.merge(df, df_new, how='inner', on=[0,1,2,3,4,5])
    fasta = read_fasta_to_dic(fasta_path)

    all_df['ref']=all_df.apply(lambda x: fasta[x[0]][x[1]] if x[5] == '+' else ref_dict[fasta[x[0]][x[1]]],axis=1)

    all_df.columns=['chrome','start','end','.','..','strand','tombo_-log10Pval',"tombo_difference",'ref']
    return all_df

def extract_xpore(xpore_path, fasta_path):
    df_original = pd.read_csv(xpore_path)
    fasta = read_fasta_to_dic(fasta_path)
    xpore_columns = df_original.columns
    for item in xpore_columns:
        if 'pval' in item:
            pval_column = item
        if 'diff_mod_rate' in item:
            diffmod_column = item
    df_original["-log10_Pval"]=df_original[pval_column].apply(lambda x: -np.log10(x))
    df_original["ref"]=df_original.apply(lambda x: x["kmer"][2], axis=1)
    df_original = df_original[['id', 'position', '-log10_Pval', diffmod_column,'ref','kmer']]
    df_original["end_position"]=df_original["position"]+1
    df_original['.']='.'
    df_original['..'] = '.'
    df_original['strand'] = df_original.apply(lambda x:'+' if x['ref'] == fasta[x['id']][x['position']] else '-',axis=1)
    df_original = df_original[['id', 'position',"end_position",'..','.','strand', '-log10_Pval', diffmod_column,'ref','kmer']]
    df_original.columns=['chrome','start','end','.','..','strand','xPore_-log10Pval',"xPore_diff_mod_rate",'ref','kmer']
    return df_original

def main(args):
    result_path = args.output
    xpore_path = args.xpore_result
    tombo_path = args.tombo_result
    fasta_path = args.ref
    if not os.path.exists(result_path):
        os.makedirs(result_path)
        print("Created the result folder")
    else:
        print("Result folder existed, will be overwrite")
    print("Start to loading tombo and xPore result . . .")
    xpore_result=extract_xpore(xpore_path,fasta_path)
    tombo_result=extract_tombo(tombo_path,fasta_path)
    print("Filtering . . .")
    final_df=pd.merge(xpore_result,tombo_result,how='outer',on=['chrome','start','end','.','..','strand','ref'])
    final_df.to_csv(result_path+"/merged_all.bed", sep='\t', index=None)
    # final_df=final_df[final_df['ref']=='A']
    result_df =final_df.copy(deep=True)
    result_df = result_df[(result_df['xPore_-log10Pval']>7) & (result_df['xPore_diff_mod_rate'].abs()>0.1)]
    result_df = result_df[(result_df['tombo_-log10Pval']>2) & (result_df['tombo_difference'].abs()>0.06)]
    result_df.to_csv(result_path+"/intersection.bed", sep='\t', index=None,header=None)
    # result_df.to_csv("/t1/zhguo/Data/sars-cov-2_result/intersection.bed", sep='\t', index=None)
    final_df['new'] = final_df.apply(lambda x:x['chrome']+"_"+ str(x['start'])+"_"+x['strand'] ,axis=1)
    xpore_pos=set(final_df[(final_df['xPore_-log10Pval']>7) & (final_df['xPore_diff_mod_rate'].abs()>0.1)]['new'].values)
    tombo_pos=set(final_df[(final_df['tombo_-log10Pval']>2) & (final_df['tombo_difference'].abs()>0.06)]['new'].values)
    print("Plotting . . .")
    my_dpi=300
    font2 = { 'size': 12}
    plt.rc('font', **font2)
    plt.figure(figsize=(4,4), dpi=my_dpi)#控制图尺寸的同时，使图高分辨率（高清）显示
    out=venn2(subsets = [xpore_pos,tombo_pos], #绘图数据集
            set_labels = ('xPore', 'tombo'), #设置组名
            set_colors=("#9790bd","#9ab4db"),
            alpha=0.8,#透明度
            normalize_to=1.0,#venn图占据figure的比例，1.0为占满
           )

    plt.title("Intersections")
    plt.savefig(result_path+"/intersection.pdf")
    print("Finished !")
    # plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tombo_result", default='sample',help="suffix of tombo result")
    parser.add_argument("--xpore_result", default="xpore_result",help="result folder from xPore")
    parser.add_argument("--ref", default='/t1/zhguo/Data/Alpha_virus_data/SINV/SINV_Toto1101.fa', help='reference path')
    parser.add_argument("--output", default='intersection_result', help='result path')
    args = parser.parse_args()
    main(args)