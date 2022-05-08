# classify reactions

import json
import pandas as pd
tf = open(r'D:\model_research\yeast_GEM_multi_omics_analysis\regulation_analysis_for_yeast9\graph\output\sample_allreaction.json','r')
sam = json.load(tf)
sample = {}
for i, j in zip(sam.keys(), sam.values()):
    if j != 'nan':
        sample[i] = j

cla = pd.read_table(r'regulation_analysis_for_yeast9/graph/data/Rxn_unique_subsystem.tsv')

cladict = dict(zip(cla.loc[:, 'ID'].values.tolist(),
                   cla.loc[:, 'subsystem_unique_@2021_12'].values.tolist()
                   )
               )

# 0.2, 0.5, inf
cla02 = {}
cla05 = {}
clainf = {}
for k, v in zip(sample.keys(), sample.values()):
    if 0.2 > v > -0.2:
        cla02[k] = v
    if -0.5 < v <= -0.2 or 0.2 <= v < 0.5:
        cla05[k] = v
    if -10 < v <= -0.5 or 0.5 <= v < 10:
        clainf[k] = v

def cla_sub(cla):
    cla_sub = {}
    for c in cla.keys():
        try:
            cla_sub[c] = cladict[c]
        except KeyError:
            continue
    count = pd.value_counts(list(cla_sub.values()))
    return count

count02 = cla_sub(cla02)
count05 = cla_sub(cla05)
countinf = cla_sub(clainf)
