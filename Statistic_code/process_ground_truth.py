import json
import pandas as pd

def ground_truth_pocessing(barcode_count,sample_name):
    # read in ground truth file and select cells whose barcodea are in barcode file
    f = open('../00_data/03_barcode/merge{}_{}_barcodes_rm-1.tsv'.format(barcode_count,sample_name))
    line = f.readline()
    d = {}
    while line:
        barcode = line.strip('\n') + '-1'
        d[barcode] = {}
        line = f.readline()
    f.close()

    df = pd.read_csv('../00_data/05_ground_truth_update/{}_ground_truth.csv'.format(sample_name))
    ground_truth_df = df.query('cell_id in @d')

    # split doublet
    ground_truth_df_1 = ground_truth_df.iloc[:,[0,1,3,5,7,9,11,13,15]]
    ground_truth_df_2 = ground_truth_df.iloc[:,[0,2,4,6,8,10,12,14,16]]
    ground_truth_df_1.columns = ['cell_id','IR_VJ_junction','IR_VDJ_junction','IR_VJ_junction_aa',
                                 'IR_VDJ_junction_aa','IR_VJ_j_call','IR_VDJ_j_call','IR_VJ_v_call',
                                 'IR_VDJ_v_call']
    ground_truth_df_2.columns = ['cell_id','IR_VJ_junction','IR_VDJ_junction','IR_VJ_junction_aa',
                                 'IR_VDJ_junction_aa','IR_VJ_j_call','IR_VDJ_j_call','IR_VJ_v_call',
                                 'IR_VDJ_v_call']
    ground_truth_df = pd.concat([ground_truth_df_1, ground_truth_df_2], ignore_index=True)
    ground_truth_df = ground_truth_df.dropna(how='all',subset=ground_truth_df.columns[1:])

    # split info from light and heavy chains
    ground_truth_df_VJ = ground_truth_df.iloc[:,[0,1,3,5,7]]
    ground_truth_df_VDJ = ground_truth_df.iloc[:,[0,2,4,6,8]]
    ground_truth_df_VJ.columns = ['cell_id','CDR3nt','CDR3aa','J','V']
    ground_truth_df_VDJ.columns = ['cell_id','CDR3nt','CDR3aa','J','V']

    ground_truth_df_VJ = ground_truth_df_VJ.dropna(how='any',subset=ground_truth_df_VJ.columns[1:])
    ground_truth_df_VDJ = ground_truth_df_VDJ.dropna(how='any',subset=ground_truth_df_VDJ.columns[1:])

    return ground_truth_df_VJ, ground_truth_df_VDJ


def query_ground_truth(sample_name, barcode_count):
    s = sample_name
    b = barcode_count
    with open('processed_ground_truth/{}_merge{}.json'.format(s,b),'r',encoding='utf-8') as f:
        ground_truth = json.load(f)
    f.close()
    ground_truth_df_VJ = pd.DataFrame.from_dict(ground_truth['VJ'])
    ground_truth_df_VDJ = pd.DataFrame.from_dict(ground_truth['VDJ'])
    return ground_truth_df_VJ, ground_truth_df_VDJ


if __name__ == '__main__':
    for sample_name in ['CPIc_C1','CPIc_C2','CPIc_C3','CPIc_C4','CPIc_C5','BC09_TUMOR1','BC09_TUMOR2','BC10_TUMOR1']:
        for barcode_count in [100,500,1000]:
            ground_truth_df_VJ, ground_truth_df_VDJ = ground_truth_pocessing(barcode_count,sample_name)
            ground_truth_df_VJ_dict = ground_truth_df_VJ.to_dict()
            ground_truth_df_VDJ_dict = ground_truth_df_VDJ.to_dict()
            ground_truth_dict = {'VJ':ground_truth_df_VJ_dict, 'VDJ':ground_truth_df_VDJ_dict}
            
            with open('processed_ground_truth/{}_merge{}.json'.format(sample_name,barcode_count), 'w') as f:
                f.write(json.dumps(ground_truth_dict))
            f.close()




