########################################
# GO enrichment visualization
# "ile_GO.txt" and "phe_GO.txt" are DAVID GO enrichment results
# get figure "GO02bar.jpg"(PCC > 0.2) and "GO-02bar.jpg"(PCC < -0.2)

import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

os.chdir('..')
# prepare figure data
goile = pd.read_table('../N_lim/output/ile_GO.txt')
gophe = pd.read_table('../N_lim/output/phe_GO.txt')
gosame = pd.read_table('../N_lim/output/Phe_Ile_same.txt')
ile4 = goile[goile['PValue'] < 0.0001]
phe4 = gophe[gophe['PValue'] < 0.0001]
same4 = gosame[gophe['PValue'] < 0.0001]
for i in ile4.index.values.tolist():
    ile4.loc[i, 'Term'] = ile4.loc[i, 'Term'].split('~')[1]
for i in phe4.index.values.tolist():
    phe4.loc[i, 'Term'] = phe4.loc[i, 'Term'].split('~')[1]
for i in same4.index.values.tolist():
    same4.loc[i, 'Term'] = same4.loc[i, 'Term'].split('~')[1]
# draw
def num2color(values, cmap):
    """Map values to colors"""
    norm = mpl.colors.Normalize(vmin=np.min(values), vmax=np.max(values))
    cmap = mpl.cm.get_cmap(cmap)
    return [cmap(norm(val)) for val in values]


colors = num2color(ile4.loc[:, 'PValue'],
                   "YlOrRd")
font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 20,
         }
font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 18,
         }
font3 = {'family': 'Arial',
         'weight': 'normal',
         'size': 16,
         }
plt.subplots(figsize=(8, 6))
plt.title('Isoleucine',
          fontdict=font1)
plt.barh(y=ile4.Term,
         color=colors,
         width=ile4.Count,
         edgecolor='#000000')
plt.xlabel('Gene number',
           fontdict=font2)
plt.yticks(fontproperties='Arial',
           size=14)
plt.xticks(fontproperties='Arial',
           size=16)
# plt.tight_layout()
plt.savefig(
    '../N_lim/output/ile_GO_bar.jpg',
    dpi=800
)
plt.show()

######################################
# fu
colors = num2color(phe4.loc[:, 'PValue'],
                   "YlOrRd")

plt.subplots(figsize=(8, 6))
plt.title('Phenylalanine',
          fontdict=font1)
plt.barh(y=phe4.Term,
         color=colors,
         width=phe4.Count,
         edgecolor='#000000')
plt.xlabel('Gene number',
           fontdict=font2)
plt.yticks(fontproperties='Arial',
           size=14)
plt.xticks(fontproperties='Arial',
           size=16)
# plt.tight_layout()
plt.savefig(
    '../N_lim/output/phe_GO_bar.jpg',
    dpi=600
)
plt.show()

#################################################
# same
plt.subplots(figsize=(8, 6))
plt.title('Shared genes',
          fontdict=font2)
plt.barh(y=same4.Term,
         color=colors,
         width=same4.Count)
plt.xlabel('Gene number',
           fontdict=font2)
plt.yticks(fontproperties='Arial',
           size=14)
plt.xticks(fontproperties='Arial',
           size=16)
plt.tight_layout()
plt.savefig(
    '../N_lim/output/same.jpg',
    dpi=800
)
plt.show()