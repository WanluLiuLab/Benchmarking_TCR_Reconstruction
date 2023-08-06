import pandas as pd 
import json
from process_ground_truth import query_ground_truth
from comparison import comparing_with_ground_truth_VJonly, comparing_with_ground_truth_CDR3only
import re

statistics_df = pd.DataFrame({'software':[],'data_type':[],'sample':[],'pseudo_count':[],'statistic':[],'value':[]})
# statistics_df_cdr3 = pd.DataFrame({'software':[],'data_type':[],'sample':[],'pseudo_count':[],'statistic':[],'value':[]})

def extract_TCR_info(filename):
    with open(filename) as file:
        content = file.readlines()

    segments = []
    i = 0
    while i < len(content):
        if "V segment" in content[i]:
            v_match = re.search(r"V segment:\t(.+)", content[i])
            v_segment = v_match.group(1)
            if "J segment" in content[i+1]:
                j_match = re.search(r"J segment:\t(.+)", content[i+1])
                j_segment = j_match.group(1)
                segments.append({"V segment": v_segment, "J segment": j_segment})
                i += 1
            elif "J segment" in content[i+2]:
                j_match = re.search(r"J segment:\t(.+)", content[i+2])
                j_segment = j_match.group(1)
                segments.append({"V segment": v_segment, "J segment": j_segment})
                i += 2
        i += 1
    file.close()
    if not segments:
        df = pd.DataFrame({'V':[],'J':[]})
    else:
        df = pd.DataFrame(segments).rename(columns={'V segment':'V','J segment':'J'})
    return df


software = 'tracer'
data_type = 'fastq'
samples = ['CPIc_C1','CPIc_C2','CPIc_C3','CPIc_C4','CPIc_C5']#,'BC09_TUMOR1','BC09_TUMOR2','BC10_TUMOR1']
counts = [100,500,1000]
lstat = ['VDJ_accuracy','VJ_accuracy','VDJ_positive_rate','VJ_positive_rate']

for s,c in [(x,y) for x in samples for y in counts]:
    print('comparing sample: '+s+' count: '+str(c)+' input type: '+data_type+'...')
    # get result
    filename = './data/05_tracer/01_PE_output/{}/{}_merge{}/unfiltered_TCR_seqs/unfiltered_TCRs.txt'.format(s,s,c)
    result = extract_TCR_info(filename)
    if not result.empty:
        for i in ['V', 'J']:
            result[i] = result[i].str.replace('\*.*', '', regex = True)  
    # get ground truth
    ground_truth_df_VJ,ground_truth_df_VDJ = query_ground_truth(s, c)
    # compare
    a_1,a_2,p_1,p_2 = comparing_with_ground_truth_VJonly(result, ground_truth_df_VJ, ground_truth_df_VDJ)
    lvalue = [a_1,a_2,p_1,p_2]
    
    for i in range(4):
        new_record = pd.Series({'software':software,
                                'data_type':data_type,
                                'sample':s,
                                'pseudo_count':c,
                                'statistic':lstat[i],
                                'value':lvalue[i]
                            })   
        statistics_df = pd.concat([statistics_df, new_record.to_frame().T], ignore_index=True)
    lstat2 = ['VDJ_accuracy','VJ_accuracy','VDJ_positive_rate','VJ_positive_rate']
    lvalue2 = list(comparing_with_ground_truth_CDR3only(result, ground_truth_df_VJ, ground_truth_df_VDJ))
    statistics_df_cdr3 = append_stat_df(software,data_type,s,c,lstat2,lvalue2,statistics_df_cdr3)


software = 'tracer'
data_type = 'fastq'
samples = ['BC09_TUMOR1','BC09_TUMOR2','BC10_TUMOR1']
counts = [100,500,1000]
lstat = ['VDJ_accuracy','VJ_accuracy','VDJ_positive_rate','VJ_positive_rate']

for s,c in [(x,y) for x in samples for y in counts]:
    print('comparing sample: '+s+' count: '+str(c)+' input type: '+data_type+'...')
    # get result
    filename = './data/05_tracer/01_SE_output/{}/{}_merge{}/unfiltered_TCR_seqs/unfiltered_TCRs.txt'.format(s,s,c)
    result = extract_TCR_info(filename)
    if not result.empty:
        for i in ['V', 'J']:
            result[i] = result[i].str.replace('\*.*', '', regex = True)  
    # get ground truth
    ground_truth_df_VJ,ground_truth_df_VDJ = query_ground_truth(s, c)
    # compare
    a_1,a_2,p_1,p_2 = comparing_with_ground_truth_VJonly(result, ground_truth_df_VJ, ground_truth_df_VDJ)
    lvalue = [a_1,a_2,p_1,p_2]
    
    for i in range(4):
        new_record = pd.Series({'software':software,
                                'data_type':data_type,
                                'sample':s,
                                'pseudo_count':c,
                                'statistic':lstat[i],
                                'value':lvalue[i]
                            })   
        statistics_df = pd.concat([statistics_df, new_record.to_frame().T], ignore_index=True)
    lstat2 = ['VDJ_accuracy','VJ_accuracy','VDJ_positive_rate','VJ_positive_rate']
    lvalue2 = list(comparing_with_ground_truth_CDR3only(result, ground_truth_df_VJ, ground_truth_df_VDJ))
    statistics_df_cdr3 = append_stat_df(software,data_type,s,c,lstat2,lvalue2,statistics_df_cdr3)


statistics_df.to_csv('/mnt/volume3/trn/04_benchmark/11_pseudo_bulkrnaseq_new/stat_result/{}_df.csv'.format(software),sep = '\t')
# statistics_df.to_csv('/mnt/volume3/trn/04_benchmark/11_pseudo_bulkrnaseq_new/stat_result/{}_df_cdr3.csv'.format(software),sep = '\t')
print('Run successfully. Result has been written to file comparison_test_df.csv in dir stat_result/{}_df.csv'.format(software))