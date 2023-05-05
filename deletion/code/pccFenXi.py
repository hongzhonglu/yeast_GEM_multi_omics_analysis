############################################
# File "Glucose_fluxDataset.csv" is from "Run.m"
# get figure "PCCfenxi.jpg" and "highmediumlow.jpg" and file "SubsystemAnalysis.xlsx"
# "SubsystemAnalysis.xlsx" is for DAVID GO enrichment

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import cobra
from scipy import stats
import re

os.chdir('..')
# prepare figure data
sub = pd.read_table('../deletion/data/uniqueSubsystems.tsv')
sub.index = sub.loc[:, 'ID']
flux = pd.read_csv('../deletion/output/Glucose_fluxDataset.csv')
flux.index = flux.loc[:, 'Var1']
flux = flux.drop(['Var1'], axis=1)
model = cobra.io.read_sbml_model('../deletion/data/yeast-GEM.xml')
trans = pd.read_csv('../deletion/data/SingleKOExp.csv')
huo_tran = pd.DataFrame(index=trans.loc[:, 'geneName'],
                        columns=flux.columns)
for t in flux.columns.values.tolist():
    huo_tran.loc[:, t] = trans.loc[:, t].values.tolist()

gene_l = []
for ge in model.genes:
    if ge.id in trans.loc[:, 'geneName'].values.tolist():
        gene_l.append(ge)

genelist = [j.id for j in gene_l]

pdata = {}
PCC = {}
for g in gene_l:
    temprxn = []
    for r in g.reactions:
        temprxn.append(r.id)
    temppcc = []
    pflux = 0
    ptran = 0
    for rx in temprxn:
        tempflux = flux.loc[rx, :]
        temptran = huo_tran.loc[g.id, :]
        r = stats.ttest_ind(tempflux, temptran)
        if np.sum(tempflux) != 0:
            pcc = np.corrcoef(tempflux, temptran)
        else:
            pcc = 0
        try:
            temppcc.append(pcc[0, 1])
        except TypeError:
            temppcc.append(pcc)
    # for genes related to more than one reations, get the PCC mean
    PCC[g.id] = np.mean(temppcc)
    pdata[g.id] = r.pvalue
PCC = pd.Series(PCC)
PCC_abs = abs(PCC)
PCC = PCC.sort_values()
PCC_abs = PCC_abs.sort_values(ascending=False)
PCC_clear = PCC.drop(PCC[(PCC[:] < 0.00001) & (PCC[:] > -0.00001)].index)
PCC.to_excel('../deletion/output/PCC_gene.xlsx')

# draw
font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 23,
         }
font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 16,
         }
font3 = {'family': 'Arial',
         'weight': 'normal',
         'size': 12,
         }
ax = plt.gca()
plt.subplots(figsize=(10, 8))
plt.bar(PCC_clear.index.values.tolist(),
        PCC_clear,
        width=1,
        color='#4E256B')
plt.xticks(color='w')
ax.axes.xaxis.set_visible(False)
plt.ylabel('PCC',
           fontdict=font2)
plt.xlabel('Gene Name',
           fontdict=font2)
plt.xticks([])
plt.yticks(fontproperties='Arial',
           size=12)
plt.savefig('../deletion/output/PCCfenxi.jpg')
plt.show()

# analysis subsystem

model = cobra.io.read_sbml_model('../deletion/data/yeast-GEM.xml')
susana = pd.DataFrame(index=PCC_abs.index,
                      columns=['PCC', 'rxn', 'subsystem'])
susana.loc[:, 'PCC'] = PCC_abs.loc[:]
modelgene = [mg.id for mg in model.genes]
for g in PCC_abs.index.values.tolist():
    temprxn = []
    tempsub = []
    if g in modelgene:
        for r in model.genes[modelgene.index(g)].reactions:
            temprxn.append(r.id)
            temprxn = list(set(temprxn))
        # susana.loc[g, 'rxn'] = '+'.join(temprxn)  for highmediumlow.jpg
        susana.loc[g, 'rxn'] = temprxn[0]  # for Stacked bar
        for rr in temprxn:
            if rr in sub.index.values.tolist():
                tempsub.append(sub.loc[rr, 'subsystem_unique_@2021_12'])
                tempsub = list(set(tempsub))
        # susana.loc[g, 'subsystem'] = '+'.join(tempsub)
        try:
            susana.loc[g, 'subsystem'] = tempsub[0]
        except IndexError:
            pass

high = susana[susana['PCC'] >= 0.4]
highcount = high.value_counts('subsystem')
medium = susana[(susana['PCC'] >= 0.2) & (susana['PCC'] < 0.4)]
mediumcount = medium.value_counts('subsystem')
low = susana[susana['PCC'] < 0.2]
lowcount = low.value_counts('subsystem')

print("high:" + str(len(high)) +
      '\nmedium:' + str(len(medium)) +
      '\nlow:' + str(len(low)))

# draw
names = ['|PCC| ≥ 0.4', '0.2 ≤ |PCC| < 0.4', '|PCC| < 0.2']
size = [len(high), len(medium), len(low)]
plt.subplots(figsize=(8, 6))
my_circle = plt.Circle((0, 0), 0.7, color='white')

plt.pie(size,
        #labels=names,
        colors=['#3491F7', '#3A51E0', '#5E34F7'],
        wedgeprops={'linewidth': 3, 'edgecolor': 'white'},
        textprops={'fontsize': 16})
p = plt.gcf()
p.gca().add_artist(my_circle)
plt.yticks(fontproperties='Arial')
plt.tight_layout()
plt.savefig('../deletion/output/highmediumlow.png',
            dpi=600,
            transparent=True)
plt.show()

# get SubsystemAnalysis.xlsx
susana_2 = pd.DataFrame(index=PCC.index,
                        columns=['PCC', 'rxn', 'subsystem', 'Standard_name'])
susana_2.loc[:, 'PCC'] = PCC.loc[:]
modelgene = [mg.id for mg in model.genes]
SGD = pd.read_table('deletion/data/SGDgeneNames.tsv')
SGD.index = SGD.loc[:, 'Systematic_name']
for g in PCC.index.values.tolist():
    temprxn = []
    tempsub = []
    if g in modelgene:
        for r in model.genes[modelgene.index(g)].reactions:
            temprxn.append(r.id)
            temprxn = list(set(temprxn))
        # susana.loc[g, 'rxn'] = '+'.join(temprxn)  for highmediumlow.jpg
        susana_2.loc[g, 'rxn'] = temprxn[0]  # for Stacked bar
        for rr in temprxn:
            if rr in sub.index.values.tolist():
                tempsub.append(sub.loc[rr, 'subsystem_unique_@2021_12'])
                tempsub = list(set(tempsub))
        # susana.loc[g, 'subsystem'] = '+'.join(tempsub)
        try:
            susana_2.loc[g, 'subsystem'] = tempsub[0]
        except IndexError:
            pass

for pr in susana_2.index.values.tolist():
    try:
        susana_2.loc[pr, 'Standard_name'] = SGD.loc[pr, 'Standard_name']
    except KeyError:
        pass

susana_2.to_excel('../deletion/output/SubsystemAnalysis.xlsx')
print('finish')
