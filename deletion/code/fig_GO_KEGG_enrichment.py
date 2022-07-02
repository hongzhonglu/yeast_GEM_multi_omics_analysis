########################################
# GO enrichment visualization
# "PCC-02_GO.txt" and "PCC-02_GO.txt" are DAVID GO enrichment results
# get figure "GO02bar.jpg"(PCC > 0.2) and "GO-02bar.jpg"(PCC < -0.2)

import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

os.chdir('..')
# prepare figure data
gofu = pd.read_table('../deletion/output/PCC-02_GO.txt')
gozheng = pd.read_table('../deletion/output/PCC02_GO.txt')
zheng5 = gozheng[gozheng['PValue'] < 0.01]
fu5 = gofu[gofu['PValue'] < 0.01]
for i in zheng5.index.values.tolist():
    zheng5.loc[i, 'Term'] = zheng5.loc[i, 'Term'].split('~')[1]
for i in fu5.index.values.tolist():
    fu5.loc[i, 'Term'] = fu5.loc[i, 'Term'].split('~')[1]

# draw
def num2color(values, cmap):
    """Map values to colors"""
    norm = mpl.colors.Normalize(vmin=np.min(values), vmax=np.max(values))
    cmap = mpl.cm.get_cmap(cmap)
    return [cmap(norm(val)) for val in values]


colors = num2color(zheng5.loc[:, 'PValue'],
                   "RdBu")
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
plt.subplots(figsize=(10, 8))
plt.barh(y=zheng5.Term,
         color=colors,
         width=zheng5.Count)
plt.xlabel('Gene number',
           fontdict=font2)
plt.yticks(fontproperties='Arial',
           size=12)
plt.yticks(fontproperties='Arial',
           size=12)
plt.tight_layout()
plt.savefig(
    '../deletion/output/GO02bar.jpg',
    dpi=600
)
plt.show()

######################################
# fu
colors = num2color(fu5.loc[:, 'PValue'],
                   "RdBu")

plt.subplots(figsize=(10, 8))
plt.barh(y=fu5.Term,
         color=colors,
         width=fu5.Count)
plt.xlabel('Gene number',
           fontdict=font2)
plt.yticks(fontproperties='Arial',
           size=12)
plt.yticks(fontproperties='Arial',
           size=12)
plt.tight_layout()
plt.savefig(
    '../deletion/output/GO-02bar.jpg',
    dpi=600
)
plt.show()
