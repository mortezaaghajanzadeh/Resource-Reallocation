#%%
import pandas as pd
data_path = r'data/raw/'
df1 = pd.read_excel(r"C:\Github\Resource-Reallocation\data\raw\count.xlsx",index_col=0)
df2 = pd.read_excel(r"C:\Github\Resource-Reallocation\data\raw\count.xlsx",sheet_name='Sheet2',index_col=0)
def find_row(l,threshold):
    for i in l:
        if i <= threshold:
            return True
    return False
df1['status'] = df1.apply(lambda x: find_row(x,4),axis=1)
df2['status'] = df2.apply(lambda x: find_row(x,4),axis=1)

df1['min'] = df1.apply(lambda x: x.iloc[:-1].min(),axis=1)
df2['min'] = df2.apply(lambda x: x.iloc[:-1].min(),axis=1)
df2 = df2[['status','min']].copy()
df1 = df1[['status','min']].copy()
#%%
df2[['status','min']].head(100)


#%%