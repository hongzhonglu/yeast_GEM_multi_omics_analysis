import cobra
import pandas as pd
from cobra.flux_analysis import double_gene_deletion

model = cobra.io.read_sbml_model(r'synthetic lethal/data/yeast-GEM.xml')
exp = pd.read_table(r'synthetic lethal/data/SGA_NxN_clustered.tsv')
exp.index = exp.loc[:, 'A']
lie = exp.index.values.tolist()
hang = exp.columns.values.tolist()
modelgene = [g.name for g in model.genes]

pair1 = [p1 for p1 in lie if p1 in modelgene]
pair2 = [p2 for p2 in hang if p2 in modelgene]

deletion = double_gene_deletion(model,
                         gene_list1=pair1,
                         gene_list2=pair2)
deletion_s = deletion.loc[:, 'ids'].values.tolist()
deletion_l = [list(d) for d in deletion_s]
simugrow = deletion.loc[:, 'growth'].values.tolist()
simu_sl = []
zeroflux = [i for i in simugrow if i < 0.000001]
num = 0
for n in deletion_l:
    if len(n) == 2:
        num += 1
print(num)

ex = []
for r1 in pair1:
    for r2 in pair2:
        if exp.loc[r1, r2] <= -0.35:
            ex.append('SL')
        elif exp.loc[r1, r2] > -0.35:
            ex.append('non_SL')
        else:
            ex.append('nan')
zeroex = 'SL' in ex
print('总计' + str(num) + '对双敲结果\n全部为不致死')
'''''
for i in deletion_l:
    ind = deletion_l.index(i)
    for r1 in pair1:
        for r2 in pair2:
            if r1 in i and r2 in i and simugrow[ind] < 0.000001:
                simu_sl.append('SL')
                if exp.loc[r1, r2] <= -0.35:
                    ex.append('SL')
                elif exp.loc[r1, r2] > -0.35:
                    ex.append('non_SL')
                else:
                    ex.append('nan')
            if r1 in i and r2 in i and simugrow[ind] > 0.000001:
                simu_sl.append('SS')
                if exp.loc[r1, r2] <= -0.35:
                    ex.append('SL')
                elif exp.loc[r1, r2] > -0.35:
                    ex.append('non_SL')
                else:
                    ex.append('nan')

tp = 0
tn = 0
fp = 0
fn = 0
tempex = []
tempsimu = []
for t in range(len(ex)):
    if ex[t] != 'nan':
        tempex.append(ex[t])
        tempsimu.append(simu[t])

for ff in range(len(tempex)):
    if tempsimu[ff] == 'non_SL' and tempex[ff] == 'non_SL':
        tn += 1
    if tempsimu[ff] == 'SL' and tempex[ff] == 'SL':
        tp += 1
    if tempsimu[ff] == 'SL' and tempex[ff] == 'non_SL':
        fp += 1
    if tempsimu[ff] == 'non_SL' and tempex[ff] == 'SL':
        fn += 1
acc = (tn + tp)/len(tempex)
'''''

