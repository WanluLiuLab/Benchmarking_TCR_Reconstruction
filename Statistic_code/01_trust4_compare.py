import pandas as pd 
import json
from process_ground_truth import query_ground_truth
from comparison import comparing_with_ground_truth, comparing_with_ground_truth_CDR3only
from append_stat_df import append_stat_df


statistics_df = pd.DataFrame({'software':[],'data_type':[],'sample':[],'pseudo_count':[],'statistic':[],'value':[]})
statistics_df_cdr3 = pd.DataFrame({'software':[],'data_type':[],'sample':[],'pseudo_count':[],'statistic':[],'value':[]})

software = 'trust4'
data_type = 'bam'
samples = ['CPIc_C1','CPIc_C2','CPIc_C3','CPIc_C4','CPIc_C5','BC09_TUMOR1','BC09_TUMOR2','BC10_TUMOR1']
counts = [100,500,1000]
lstat = ['VDJ_accuracy','VJ_accuracy','VDJ_positive_rate','VJ_positive_rate']

for s,c in [(x,y) for x in samples for y in counts]:
    print('comparing sample: '+s+' count: '+str(c)+'input type: '+data_type+'...')
    # get result
    result = pd.read_csv('./data/01_trust4/01_bam_output/{}/{}_merge{}_report.tsv'.format(s,s,c),sep = '\t')
    for i in ['V', 'D', 'J', 'C']:
        result[i] = result[i].str.replace('\*.*', '', regex = True)    
    # get ground truth
    ground_truth_df_VJ,ground_truth_df_VDJ = query_ground_truth(s, c)
    a_1,a_2,p_1,p_2 = comparing_with_ground_truth(result, ground_truth_df_VJ, ground_truth_df_VDJ)
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




software = 'trust4'
data_type = 'fastq'
samples = ['CPIc_C1','CPIc_C2','CPIc_C3','CPIc_C4','CPIc_C5']#,'BC09_TUMOR1','BC09_TUMOR2','BC10_TUMOR1']
counts = [100,500,1000]
lstat = ['VDJ_accuracy','VJ_accuracy','VDJ_positive_rate','VJ_positive_rate']

for s,c in [(x,y) for x in samples for y in counts]:
    print('comparing sample: '+s+' count: '+str(c)+'input type: '+data_type+'...')
    # get result
    result = pd.read_csv('./data/01_trust4/01_PE_output/{}/{}_merge{}_report.tsv'.format(s,s,c),sep = '\t')
    for i in ['V', 'D', 'J', 'C']:
        result[i] = result[i].str.replace('\*.*', '', regex = True)    
    # get ground truth
    ground_truth_df_VJ,ground_truth_df_VDJ = query_ground_truth(s, c)
    a_1,a_2,p_1,p_2 = comparing_with_ground_truth(result, ground_truth_df_VJ, ground_truth_df_VDJ)
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

software = 'trust4'
data_type = 'fastq'
samples = ['BC09_TUMOR1','BC09_TUMOR2','BC10_TUMOR1']
counts = [100,500,1000]
lstat = ['VDJ_accuracy','VJ_accuracy','VDJ_positive_rate','VJ_positive_rate']

for s,c in [(x,y) for x in samples for y in counts]:
    print('comparing sample: '+s+' count: '+str(c)+'input type: '+data_type+'...')
    # get result
    result = pd.read_csv('./data/01_trust4/01_SE_output/{}/{}_merge{}_report.tsv'.format(s,s,c),sep = '\t')
    for i in ['V', 'D', 'J', 'C']:
        result[i] = result[i].str.replace('\*.*', '', regex = True)    
    # get ground truth
    ground_truth_df_VJ,ground_truth_df_VDJ = query_ground_truth(s, c)
    a_1,a_2,p_1,p_2 = comparing_with_ground_truth(result, ground_truth_df_VJ, ground_truth_df_VDJ)
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
statistics_df_cdr3.to_csv('/mnt/volume3/trn/04_benchmark/11_pseudo_bulkrnaseq_new/stat_result/{}_df_cdr3.csv'.format(software),sep = '\t')

print('run successfully, result has been written to file comparison_test_df.csv in dir stat_result/{}_df.csv'.format(software))





