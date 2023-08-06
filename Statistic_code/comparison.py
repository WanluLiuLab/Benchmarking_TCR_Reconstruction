import pandas as pd 
import json
from process_ground_truth import query_ground_truth


def comparing_with_ground_truth(result, ground_truth_df_VJ, ground_truth_df_VDJ, CDR3_form='nt'):
    '''
    This function use compare the difference of VDJ gene and CDR3 seq between the ground truth generate by
    TCR sequencing and pseudo bulk analysis, and calculate the true positive rate and accuracy of pesudo 
    bulk RNA seq TCR recognizing result.
    '''
 # ========================================================================================================= #    
    # These code should be altered when process result from different software
    # notice the type of result data (fastq/bam)
    
    # result = pd.read_csv('{}/{}_merge{}_report.tsv'.format(sample_name,sample_name,barcode_count),sep = '\t')
    # for i in ['V', 'D', 'J', 'C']:
    #     result[i] = result[i].str.replace('\*.*', '', regex = True)
 # ========================================================================================================= #    
    # Here, this function should take in result - (dictionary contain two Pandas dataframe with column CDR3nt,
    # V and J) directly as its parameter.
    
    # part1 true positive rate
    result = result.reset_index()
    result_cnt = len(result)
    if result_cnt == 0:
        print('No detected TCR result in this file!')
        return 0,0,0,0
    res_match = 0
    VJ_res_match = 0
    VDJ_res_match = 0
    VJ_res_total = 0
    VDJ_res_total = 0
    
    # Some software only provide CDR3 as animo acid seq. Here by changing the column name we compare the aa seq
    # using the same code as process comparing CDR3nt
    if CDR3_form == 'aa':
        result = result.rename(columns={'CDR3aa':'CDR3nt'})
        ground_truth_df_VJ = ground_truth_df_VJ.rename(columns={'CDR3aa':'CDR3nt','CDR3nt':'CDR3nt_real'})
        ground_truth_df_VDJ = ground_truth_df_VDJ.rename(columns={'CDR3aa':'CDR3nt','CDR3nt':'CDR3nt_real'})

    for i in range(result_cnt):
        # information of CDR3 sequence, V gene and J gene should be extracted
        CDR3nt = result.at[i,'CDR3nt']
        V_gene = result.at[i,'V']
        J_gene = result.at[i,'J']
        if len(V_gene) < 5:
            V_gene = 'None'
        if len(J_gene) < 5:
            J_gene = 'None'
        if V_gene[2] == 'A' or J_gene[2] == 'A': # query in VJ df
            VJ_res_total += 1
            query_df = ground_truth_df_VJ.query('CDR3nt == @CDR3nt and V == @V_gene and J == @J_gene')
            if not query_df.empty:
                res_match += 1
                VJ_res_match += 1   
        elif V_gene[2] == 'B' or J_gene[2] == 'B': # query in VJ df
            VDJ_res_total += 1
            query_df = ground_truth_df_VDJ.query('CDR3nt == @CDR3nt and V == @V_gene and J == @J_gene')
            if not query_df.empty:
                res_match += 1
                VDJ_res_match += 1
                
        # VJ_query_df = ground_truth_df_VJ.query('CDR3nt == @CDR3nt').reset_index()
        # VDJ_query_df = ground_truth_df_VDJ.query('CDR3nt == @CDR3nt').reset_index()
            

        # # The default here is that it is impossible to find a match in both dataframe
        # if not VJ_query_df.empty:
        #     if (VJ_query_df.at[0,'J']==J_gene) and (VJ_query_df.at[0,'V']==V_gene):
        #         res_match += 1
        # elif not VDJ_query_df.empty:
        #     if (VDJ_query_df.at[0,'J']==J_gene) and (VDJ_query_df.at[0,'V']==V_gene):
        #         res_match += 1
        # else:
        #     continue

    positive_rate = res_match/result_cnt
    VDJ_positive_rate = VDJ_res_match/VDJ_res_total if VDJ_res_total != 0 else 0
    VJ_positive_rate = VJ_res_match/VJ_res_total if VJ_res_total != 0 else 0
    
    # part2 accuracy
    total_truth_match = 0
    VJ_truth_match = 0
    VDJ_truth_match = 0
    
    ground_truth_df_VJ = ground_truth_df_VJ[ground_truth_df_VJ['CDR3nt'].notna()].reset_index()
    ground_truth_df_VDJ = ground_truth_df_VDJ[ground_truth_df_VDJ['CDR3nt'].notna()].reset_index()
    ground_truth_count = len(ground_truth_df_VJ) + len(ground_truth_df_VDJ)

    for i in range(len(ground_truth_df_VJ)):
        CDR3nt = ground_truth_df_VJ.at[i,'CDR3nt']
        V_gene = ground_truth_df_VJ.at[i,'V']
        J_gene = ground_truth_df_VJ.at[i,'J']
        query_df = result.query('CDR3nt == @CDR3nt and V == @V_gene and J == @J_gene')
        if not query_df.empty:
            total_truth_match += 1
            VJ_truth_match += 1

    for i in range(len(ground_truth_df_VDJ)):
        CDR3nt = ground_truth_df_VDJ.at[i,'CDR3nt']
        V_gene = ground_truth_df_VDJ.at[i,'V']
        J_gene = ground_truth_df_VDJ.at[i,'J']
        query_df = result.query('CDR3nt == @CDR3nt and V == @V_gene and J == @J_gene')
        if not query_df.empty:
            total_truth_match += 1
            VDJ_truth_match += 1

    accuracy = total_truth_match/ground_truth_count
    # print('VJ: {}/{}'.format(VJ_truth_match,len(ground_truth_df_VJ)))
    # print('VDJ: {}/{}'.format(VDJ_truth_match,len(ground_truth_df_VDJ)))
    VJ_accuracy = VJ_truth_match/len(ground_truth_df_VJ)
    VDJ_accuracy = VDJ_truth_match/len(ground_truth_df_VDJ)
    
    return VDJ_accuracy,VJ_accuracy,VDJ_positive_rate,VJ_positive_rate
    

def comparing_with_ground_truth_VJonly(result, ground_truth_df_VJ, ground_truth_df_VDJ):
    '''
    Simplified version of func comparing_with_ground_truth,
    as the software TraCeR didn't provide CDR3 info, here only V gene and J gene will be compared.
    '''
    result = result.reset_index()
    result_cnt = len(result)
    if result_cnt == 0:
        print('No detected TCR result in this file!')
        return 0,0,0,0
    res_match = 0
    
    # part1 true positive
    for i in range(result_cnt):
        V_gene = result.at[i,'V']
        J_gene = result.at[i,'J']
        if V_gene[2] == 'A': # query in VJ df
            query_df = ground_truth_df_VJ.query('V == @V_gene and J == @J_gene')
            if not query_df.empty:
                res_match += 1
        elif V_gene[2] == 'B': # query in VJ df
            query_df = ground_truth_df_VDJ.query('V == @V_gene and J == @J_gene')
            if not query_df.empty:
                res_match += 1

    positive_rate = res_match/result_cnt
    
    # part2 accuracy
    total_truth_match = 0
    VJ_truth_match = 0
    VDJ_truth_match = 0
    
    ground_truth_df_VJ = ground_truth_df_VJ[ground_truth_df_VJ['CDR3nt'].notna()].reset_index()
    ground_truth_df_VDJ = ground_truth_df_VDJ[ground_truth_df_VDJ['CDR3nt'].notna()].reset_index()
    ground_truth_count = len(ground_truth_df_VJ) + len(ground_truth_df_VDJ)

    for i in range(len(ground_truth_df_VJ)):
        V_gene = ground_truth_df_VJ.at[i,'V']
        J_gene = ground_truth_df_VJ.at[i,'J']
        query_df = result.query('V == @V_gene and J == @J_gene')
        if not query_df.empty:
            total_truth_match += 1
            VJ_truth_match += 1

    for i in range(len(ground_truth_df_VDJ)):
        V_gene = ground_truth_df_VDJ.at[i,'V']
        J_gene = ground_truth_df_VDJ.at[i,'J']
        query_df = result.query('V == @V_gene and J == @J_gene')
        if not query_df.empty:
            total_truth_match += 1
            VDJ_truth_match += 1

    accuracy = total_truth_match/ground_truth_count
    VJ_accuracy = VJ_truth_match/len(ground_truth_df_VJ)
    VDJ_accuracy = VDJ_truth_match/len(ground_truth_df_VDJ)
    
    return accuracy,VDJ_accuracy,VJ_accuracy,positive_rate


def comparing_with_ground_truth_CDR3only(result, ground_truth_df_VJ, ground_truth_df_VDJ,CDR3_form='nt'):
    '''
    Simplified version of func comparing_with_ground_truth, only compare CDR3 infomation
    '''
    result = result.reset_index()
    result_cnt = len(result)
    if result_cnt == 0:
        print('No detected TCR result in this file!')
        return 0,0,0,0

    res_match = 0
    VJ_res_match = 0
    VDJ_res_match = 0
    VJ_res_total = 0
    VDJ_res_total = 0
    
    if CDR3_form == 'aa':
        result = result.rename(columns={'CDR3aa':'CDR3nt'})
        ground_truth_df_VJ = ground_truth_df_VJ.rename(columns={'CDR3aa':'CDR3nt','CDR3nt':'CDR3nt_real'})
        ground_truth_df_VDJ = ground_truth_df_VDJ.rename(columns={'CDR3aa':'CDR3nt','CDR3nt':'CDR3nt_real'})

    # part1 true positive
    for i in range(result_cnt):
        CDR3nt = result.at[i,'CDR3nt']
        V_gene = result.at[i,'V']
        J_gene = result.at[i,'J']
        if len(V_gene) < 5:
            V_gene = 'None'
        if len(J_gene) < 5:
            J_gene = 'None'
        if V_gene[2] == 'A' or J_gene[2] == 'A': # query in VJ df
            VJ_res_total += 1
            query_df = ground_truth_df_VJ.query('CDR3nt == @CDR3nt')
            if not query_df.empty:
                res_match += 1
                VJ_res_match += 1   
        elif V_gene[2] == 'B' or J_gene[2] == 'B': # query in VJ df
            VDJ_res_total += 1
            query_df = ground_truth_df_VDJ.query('CDR3nt == @CDR3nt')
            if not query_df.empty:
                res_match += 1
                VDJ_res_match += 1

    positive_rate = res_match/result_cnt
    # if VDJ_res_total != 0:
    #     VDJ_positive_rate = VDJ_res_match/VDJ_res_total
    # else:
    #     VDJ_positive_rate = 0
    VDJ_positive_rate = VDJ_res_match/VDJ_res_total if VDJ_res_total != 0 else 0
    VJ_positive_rate = VJ_res_match/VJ_res_total if VJ_res_total != 0 else 0
    
    # part2 accuracy
    VJ_truth_match = 0
    VDJ_truth_match = 0
    
    ground_truth_df_VJ = ground_truth_df_VJ[ground_truth_df_VJ['CDR3nt'].notna()].reset_index()
    ground_truth_df_VDJ = ground_truth_df_VDJ[ground_truth_df_VDJ['CDR3nt'].notna()].reset_index()

    for i in range(len(ground_truth_df_VJ)):
        CDR3nt = ground_truth_df_VJ.at[i,'CDR3nt']
        query_df = result.query('CDR3nt == @CDR3nt')
        if not query_df.empty:
            VJ_truth_match += 1

    for i in range(len(ground_truth_df_VDJ)):
        CDR3nt = ground_truth_df_VDJ.at[i,'CDR3nt']
        query_df = result.query('CDR3nt == @CDR3nt')
        if not query_df.empty:
            VDJ_truth_match += 1

    VJ_accuracy = VJ_truth_match/len(ground_truth_df_VJ)
    VDJ_accuracy = VDJ_truth_match/len(ground_truth_df_VDJ)
    
    return VDJ_accuracy,VJ_accuracy,VDJ_positive_rate,VJ_positive_rate


def testing():
    print("==="*10+'testing...'+"==="*10)
    # directly run this script for testing purpose.
    statistics_df = pd.DataFrame({'software':[],'data_type':[],'sample':[],'pseudo_count':[],
                                  'statistic':[],'value':[]})

    # 01_trust4
    software = 'trust4'
    data_type = 'bam'
    samples = ['CPIc_C1','CPIc_C2','CPIc_C3','CPIc_C4','CPIc_C5']
    counts = [100,500,1000]
    lstat = ['accuracy','VDJ_accuracy','VJ_accuracy','positive_rate']

    for s,c in [(x,y) for x in samples for y in counts]:
        print('comparing '+s+' and '+str(c)+'...')
        # s:sample_name; c:barcode_count
        # get result
        result = pd.read_csv('./data/01_trust4/01_bam_output/{}/{}_merge{}_report.tsv'.format(s,s,c),sep = '\t')
        for i in ['V', 'D', 'J', 'C']:
            result[i] = result[i].str.replace('\*.*', '', regex = True)    
        
        # get ground truth
        ground_truth_df_VJ,ground_truth_df_VDJ = query_ground_truth(s, c)
        # a:accuracy; a_1:VDJ_accuracy; a_2:VJ_accuracy; p:positive_rate
        a,a_1,a_2,p = comparing_with_ground_truth(result, ground_truth_df_VJ, ground_truth_df_VDJ)
        lvalue = [a,a_1,a_2,p]
        
        for i in range(4):
            new_record = pd.Series({'software':software,
                                    'data_type':data_type,
                                    'sample':s,
                                    'pseudo_count':c,
                                    'statistic':lstat[i],
                                    'value':lvalue[i]
                                })
            
            statistics_df = pd.concat([statistics_df, new_record.to_frame().T], ignore_index=True)

    statistics_df.to_csv('./tmp/comparison_test_df.csv',sep = '\t')
    print("==="*10+'test finish'+"==="*10)
    print('run successfully, result has been written to file comparison_test_df.csv in dir ./tmp/')


if __name__ == '__main__':
    testing()