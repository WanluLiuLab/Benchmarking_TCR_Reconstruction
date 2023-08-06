import pandas as pd

def append_stat_df(software,data_type,s,c,lstat,lvalue,df):
    assert len(lstat) == len(lvalue), 'length of lstat and l value should be equal'
    for i in range(len(lstat)):
        new_record = pd.Series({'software':software,
                                'data_type':data_type,
                                'sample':s,
                                'pseudo_count':c,
                                'statistic':lstat[i],
                                'value':lvalue[i]
                            })   
        df = pd.concat([df, new_record.to_frame().T], ignore_index=True)
    return df

