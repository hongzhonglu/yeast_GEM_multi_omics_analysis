#############################################
# predict synthetic lethal
# Data:
# A. H. Tong et al. , Science 294, 2364 (2001)
# A. S. Goehring et al. , Mol. Biol. Cell , 14,1501 (2003)
# K. Kozminski et al. , Mol. Biol. Cell , 14,4958 (2003)
# N. Krogan et al. , Mol. Cell Biol. , 23, 4207 (2003)
# M. Bellaoui et al. , EMBO , 22, 4304 (2003)
# D. Huang et al. , Mol. Cell Biol. , 22, 5076 (2003)
# A. H. Tong et al. , Science (2004)


import cobra
import pandas as pd
from cobra.flux_analysis import double_gene_deletion
import os

os.chdir('..')

model = cobra.io.read_sbml_model('../synthetic lethal/data/yeast-GEM.xml')
SGD = pd.read_table('../synthetic lethal/data/SGDgeneNames.tsv')
exp = pd.read_csv('../synthetic lethal/data/science2004SI.csv')
stan_name = SGD.loc[:, 'Standard_name'].values.tolist()
pair2temp = exp.loc[:, 'GeneticInteraction-SystematicName'].values.tolist()
pair1temp = exp.loc[:, 'QueryGene'].values.tolist()
S_or_L = exp.loc[:, 'synthetic lethal'].values.tolist()
pair1 = []
pair2 = []
modelgenelist = [g.name for g in model.genes]
model_sys = []
for s in modelgenelist:
    try:
        tmp = stan_name.index(s)
        model_sys.append(SGD.loc[tmp, 'Systematic_name'])
    except ValueError:
        model_sys.append(s)

pair1tempsys = []
for p in pair1temp:
    if p in stan_name:
        pair1tempsys.append(SGD.loc[stan_name.index(p), 'Systematic_name'])
    else:
        pair1tempsys.append('nan')
sl = []
n = 0
for p1, p2 in zip(pair1tempsys, pair2temp):
    if p1 in model_sys and p2 in model_sys:
        pair1.append(p1)
        pair2.append(p2)
        sl.append(S_or_L[n])
    n += 1

deletion = double_gene_deletion(model,
                         gene_list1=pair1,
                         gene_list2=pair2)
deletion_s = deletion.loc[:, 'ids'].values.tolist()
deletion_l = [list(d) for d in deletion_s]
simugrow = deletion.loc[:, 'growth'].values.tolist()
simu_sl = []
for r1, r2 in zip(pair1, pair2):
    for i in deletion_l:
        ind = deletion_l.index(i)
        if r1 in i and r2 in i and simugrow[ind] < 0.000001:
            simu_sl.append('SL')
        if r1 in i and r2 in i and simugrow[ind] > 0.000001:
            simu_sl.append('SS')

tn = 0
tp = 0
fp = 0
fn = 0
fp_l = []
fn_l = []
for ff in range(len(sl)):
    if simu_sl[ff] == 'SS' and sl[ff] == 'SS':
        tn += 1
    if simu_sl[ff] == 'SL' and sl[ff] == 'SL':
        tp += 1
    if simu_sl[ff] == 'SL' and sl[ff] == 'SS':
        fp += 1
        fp_l.append(pair1[ff] + ' and ' + pair2[ff])
    if simu_sl[ff] == 'SS' and sl[ff] == 'SL':
        fn += 1
        fn_l.append(pair1[ff] + ' and ' + pair2[ff])

acc = (tn + tp)/48

print(
    'accuracy:' + str(acc) +
    '\ntp:' + str(tp) +
    '\ntn:' + str(tn) +
    '\nfp:' + str(fp) +
    '\nfn:' + str(fn)
)
