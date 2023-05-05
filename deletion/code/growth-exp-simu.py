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
flux1 = pd.read_csv('../deletion/output/Glucose_fluxDataset.csv')
exp1 = pd.read_csv('../deletion/data/SingleGrowthExp.csv')
exp1.index = exp1.loc[:, 'commonName']
flux1.index = flux1.loc[:, 'Var1']
flux1 = flux1.drop(['Var1'], axis=1)
simutemp1 = flux1.loc['r_2111', :].values.tolist()
simu1 = pd.Series(simutemp1, index=flux1.columns.values.tolist())
simu_exp1 = pd.DataFrame(index=flux1.columns.values.tolist(),
                        columns=['exp', 'simu'])
for g in flux1.columns.values.tolist():
    simu_exp1.loc[g, 'simu'] = simu1.loc[g]
    simu_exp1.loc[g, 'exp'] = exp1.loc[g, 'log2relT']

# double
flux2 = pd.read_csv('../deletion/output/Doubleglucose_fluxDataset.csv')
exp2 = pd.read_csv('../deletion/data/DoubleGrowthExp.csv')
exp2.index = exp2.loc[:, 'commonName']
flux2.index = flux2.loc[:, 'Var1']
flux2 = flux2.T
simutemp2 = flux2.loc[:, 'r_2111'].values.tolist()
simu2 = pd.Series(simutemp2[1:], index=flux2.index.values.tolist()[1:])
simu_exp2 = pd.DataFrame(index=flux2.index.values.tolist()[1:],
                        columns=['exp', 'simu'])

for g in flux2.index.values.tolist()[1:]:
    simu_exp2.loc[g, 'simu'] = simu2.loc[g]
    simu_exp2.loc[g, 'exp'] = exp2.loc[g, 'log2relT']

simu_exp = pd.concat([simu_exp1, simu_exp2])

# draw
font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 23,
         }
font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 18,
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
           size=16)
plt.yticks(fontproperties='Arial',
           size=16)
plt.tight_layout()
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
