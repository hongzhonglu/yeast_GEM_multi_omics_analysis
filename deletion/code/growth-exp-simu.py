############################################
# File "Glucose_fluxDataset.csv" is from "Run.m"
# get figure growth_exp_simu.jpg
# compute PCC and p value between experiment data and model simulation


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import stats

os.chdir('..')

# prepare figure data
flux = pd.read_csv('../deletion/output/Glucose_fluxDataset.csv')
exp = pd.read_csv('../deletion/data/SingleGrowthExp.csv')
exp.index = exp.loc[:, 'commonName']
flux.index = flux.loc[:, 'Var1']
flux = flux.drop(['Var1'], axis=1)
simutemp = flux.loc['r_2111', :].values.tolist()
simu = pd.Series(simutemp, index=flux.columns.values.tolist())
simu_exp = pd.DataFrame(index=flux.columns.values.tolist(),
                        columns=['exp', 'simu'])

for g in flux.columns.values.tolist():
    simu_exp.loc[g, 'simu'] = simu.loc[g]
    simu_exp.loc[g, 'exp'] = exp.loc[g, 'log2relT']

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

plt.scatter(simu_exp.loc[:, 'simu'],
            simu_exp['exp'],
            alpha=0.3)
plt.ylabel('log2 of doubling time fold change',
           fontdict=font2)
plt.xlabel('Growth rate(h-1)',
           fontdict=font2)
plt.xticks(fontproperties='Arial',
           size=12)
plt.yticks(fontproperties='Arial',
           size=12)
plt.savefig('../deletion/output/growth_exp_simu.jpg')
plt.show()

# compute PCC and p values
pccs = np.corrcoef(np.array(simu_exp.loc[:, 'simu'].values.tolist()),
                   np.array(simu_exp['exp'].values.tolist()))
pairp = stats.ttest_rel(np.array(simu_exp.loc[:, 'simu'].values.tolist()),
                        np.array(simu_exp['exp'].values.tolist())
                        , axis=0, nan_policy='propagate')

print(
    'PCC:' + str(pccs[0, 1]) +
    '\np value:' + str(pairp.pvalue)
)
