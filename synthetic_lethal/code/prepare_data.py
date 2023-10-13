import pandas as pd
import cobra


model = cobra.io.read_sbml_model('../data/yeast-GEM.xml')
data2016 = pd.read_excel('../data/SGA_NxN_clustered.xlsx')
data2016.index = data2016.loc[:, 'A']
ind = data2016.index.values.tolist()
col = data2016.columns.values.tolist()
in_ind = [l for l in ind if l in model.genes]
in_col = [h for h in col if h in model.genes]
data2016 = data2016.fillna(-1000)


data2016_tmp1 = pd.DataFrame(columns=in_col, index=ind)
n = 0
for c in in_col:
    print(n)
    data2016_tmp1.loc[:, c] = data2016.loc[:, c]
    n += 1

m = 0
data2016_in = pd.DataFrame(columns=in_col, index=in_ind)
for i in in_ind:
    print(m)
    data2016_in.loc[i, :] = data2016.loc[i, :]
    m += 1

chongfu = [g for g in in_ind if g in in_col]
print('{} is repeated'.format(len(chongfu)))

data2016_in.to_excel('../data/data2016_in_model.xlsx')
pd.Series(in_ind).to_csv('../genelist_row.csv')
pd.Series(in_col).to_csv('../genelist_col.csv')

