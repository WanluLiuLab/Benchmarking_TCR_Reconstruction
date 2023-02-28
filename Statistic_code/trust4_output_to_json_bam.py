import os
import pandas as pd
import numpy as np
from collections import Counter
import demjson
import json
import sys
file_name=sys.argv[1] #
input_filedir=sys.argv[2] #02_bam_statistic
output_filedir=sys.argv[3] #03_bam_stat_json
cdr3_filedir=sys.argv[4] #01_bam_output
file_list=pd.read_table(file_name,header=None).iloc[:,0].values
#print(file_list)
def translate(seq):
     
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    #print(seq)
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon in table.keys():
                protein+= table[codon]
            else:
                protein+= '_'
    elif len(seq)%3 == 1:
        for i in range(0, len(seq)-1, 3):
            codon = seq[i:i + 3]
            if codon in table.keys():
                protein+= table[codon]
            else:
                protein+= '_'
    elif len(seq)%3 == 2:
        for i in range(0, len(seq)-2, 3):
            codon = seq[i:i + 3]
            if codon in table.keys():
                protein+= table[codon]
            else:
                protein+= '_'
    return protein

for filehead in file_list:
    print(filehead)
    candidate_reads_df=pd.read_table(
        "./"+input_filedir+"/"+filehead+"/01_candidatereads_stat_merge.txt",header=None)
    candidate_reads_stat_df=pd.DataFrame({
        "cell_barcode":candidate_reads_df.iloc[:,0].values,
        "candidate_reads_number":candidate_reads_df.iloc[:,1]//4,
        "total_reads_number":candidate_reads_df.iloc[:,2].values
    })
    assembled_reads_df=pd.read_table(
        "./"+input_filedir+"/"+filehead+"/02_assembledreads_stat_merge.txt",header=None)
    assembled_reads_stat_df=pd.DataFrame({
        "cell_barcode":assembled_reads_df.iloc[:,0].values,
        "assembled_reads_number":assembled_reads_df.iloc[:,1]//2
    })
    result1_df=pd.merge(candidate_reads_stat_df,assembled_reads_stat_df,
                        how='outer',on='cell_barcode')
    result1_df.to_csv(
        "./"+output_filedir+"/01_candidatereads/"+filehead+"_candidatereads.csv",index=None)
    contigreads_df=pd.read_table(
        "./"+input_filedir+"/"+filehead+"/03_contigreads_merge.txt",header=None)
    contigreads_annot_df=pd.read_table(
        "./"+input_filedir+"/"+filehead+"/04_contigreads_annot.txt",header=None,delimiter=' ')
    contigreads_annot_split_df=pd.DataFrame(
        {
            'cell_barcode':contigreads_annot_df.iloc[:,1],
            'contig_name':contigreads_annot_df.iloc[:,3],
            'contig_length':contigreads_annot_df.iloc[:,4],
            'contig_reads':contigreads_df.iloc[:,1]
        })
    contig_info_df=contigreads_annot_split_df
    contig_info_df.to_csv(
        "./"+output_filedir+"/02_contigs/"+filehead+"_contigs.csv",index=None)
    cellbarcode_ls=result1_df['cell_barcode'].values
    results_json=[]
    for cellbarcode_name in cellbarcode_ls:
        tpm_2_df=result1_df[result1_df['cell_barcode']==cellbarcode_name]
        if cellbarcode_name in contig_info_df["cell_barcode"].values:
            tpm_df=contig_info_df[contig_info_df["cell_barcode"]==cellbarcode_name]
            dict_tmp={
                'cellbarcode':cellbarcode_name,
                'total_reads_number':int(tpm_2_df['total_reads_number'].values[0]),
                'candidate_reads_number':int(tpm_2_df['candidate_reads_number'].values[0]),
                'assembled_reads_number':tpm_2_df['assembled_reads_number'].values[0],                
                'contig_number':tpm_df.shape[0],
                'contig_length':list(map(int,tpm_df['contig_length'].values)),
                'contig_reads':list(tpm_df['contig_reads'].values)
            }
            cdr3_file_path="./"+cdr3_filedir+"/"+filehead+"/"+filehead+"."+cellbarcode_name+".demultiplexed_cdr3.out"
            if os.path.getsize(cdr3_file_path)==0:
                dict_tmp['TRA_J_gene']=None
                dict_tmp['TRA_J_gene_number']=0
                dict_tmp['TRA_V_gene']=None
                dict_tmp['TRA_V_gene_number']=0
                dict_tmp['TRA_CDR3']=None
                dict_tmp['TRA_CDR3_AA']=None
                dict_tmp['TRA_CDR3_number']=0
                dict_tmp['TRB_J_gene']=None
                dict_tmp['TRB_J_gene_number']=0
                dict_tmp['TRB_V_gene']=None
                dict_tmp['TRB_V_gene_number']=0
                dict_tmp['TRB_CDR3']=None
                dict_tmp['TRB_CDR3_AA']=None
                dict_tmp['TRB_CDR3_number']=0
            else:
                cdr3_df=pd.read_table(cdr3_file_path,header=None)
                cdr3_df['chaintype_isTRB']=list(map(lambda x: "TRB" in ''.join(cdr3_df.iloc[x,2:5]),range(0,cdr3_df.shape[0])))
                cdr3_df['chaintype_isTRA']=list(map(lambda x: "TRA" in ''.join(cdr3_df.iloc[x,2:5]),range(0,cdr3_df.shape[0])))
                cdr3_TRB_df=cdr3_df[cdr3_df.chaintype_isTRB==True]
                cdr3_TRA_df=cdr3_df[cdr3_df.chaintype_isTRA==True]
                dict_tmp['TRA_J_gene']=list(cdr3_TRA_df.iloc[:,4])
                dict_tmp['TRA_J_gene_number']=len(list(cdr3_TRA_df.iloc[:,4]))
                dict_tmp['TRA_V_gene']=list(cdr3_TRA_df.iloc[:,2])
                dict_tmp['TRA_V_gene_number']=len(list(cdr3_TRA_df.iloc[:,2]))
                dict_tmp['TRA_CDR3']=list(cdr3_TRA_df.iloc[:,8])
                dict_tmp['TRA_CDR3_AA']=list(map(lambda x: translate(x),list(cdr3_TRA_df.iloc[:,8])))
                dict_tmp['TRA_CDR3_number']=len(list(cdr3_TRA_df.iloc[:,8]))
                dict_tmp['TRB_J_gene']=list(cdr3_TRB_df.iloc[:,4])
                dict_tmp['TRB_J_gene_number']=len(list(cdr3_TRB_df.iloc[:,4]))
                dict_tmp['TRB_V_gene']=list(cdr3_TRB_df.iloc[:,2])
                dict_tmp['TRB_V_gene_number']=len(list(cdr3_TRB_df.iloc[:,2]))
                dict_tmp['TRB_CDR3']=list(cdr3_TRB_df.iloc[:,8])
                dict_tmp['TRB_CDR3_AA']=list(map(lambda x: translate(x),list(cdr3_TRB_df.iloc[:,8])))
                dict_tmp['TRB_CDR3_number']=len(list(cdr3_TRB_df.iloc[:,8]))                
        else:
            dict_tmp={
                'cellbarcode':cellbarcode_name,
                'total_reads_number':int(tpm_2_df['total_reads_number'].values[0]),
                'candidate_reads_number':int(tpm_2_df['candidate_reads_number'].values[0]),
                'assembled_reads_number':tpm_2_df['assembled_reads_number'].values[0], 
                'contig_number':0,
                'contig_length':0,
                'contig_reads':None
            }
            dict_tmp['TRA_J_gene']=None
            dict_tmp['TRA_J_gene_number']=0
            dict_tmp['TRA_V_gene']=None
            dict_tmp['TRA_V_gene_number']=0
            dict_tmp['TRA_CDR3']=None
            dict_tmp['TRA_CDR3_AA']=None
            dict_tmp['TRA_CDR3_number']=0
            dict_tmp['TRB_J_gene']=None
            dict_tmp['TRB_J_gene_number']=0
            dict_tmp['TRB_V_gene']=None
            dict_tmp['TRB_V_gene_number']=0
            dict_tmp['TRB_CDR3']=None
            dict_tmp['TRB_CDR3_AA']=None
            dict_tmp['TRB_CDR3_number']=0
        results_json.append(dict_tmp)
    f=open("./"+output_filedir+"/"+filehead+"trust4_bam_results.json",'w')
    f.write(json.dumps(results_json))
    f.close()    


