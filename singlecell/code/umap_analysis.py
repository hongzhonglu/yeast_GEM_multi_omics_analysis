############################################
# apply umap to visualize classification
# get "umap.jpg"
# first run "get_pfba_flux.py" to compute flux data "pfba_flux.xlsx"

import umap
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pandas as pd
import os

os.chdir('..')

umap_data = pd.read_excel('../singlecell/output/pfba_flux.xlsx')
umap_data = umap_data.drop(0, axis=1)
umap_data.loc[:, 'metnum/20'] = umap_data.loc[:, 'metnum/20'] * 20
umap_data.loc[:, 'rxnnum/20'] = umap_data.loc[:, 'rxnnum/20'] * 20

reducer = umap.UMAP(n_components=3)
scaled_umap_data = StandardScaler().fit_transform(umap_data)
embedding = reducer.fit_transform(scaled_umap_data)


lw = 2
y = [0]*80 + [1]*83
font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 23,
         }
font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 20,
         }
font3 = {'family': 'Arial',
         'weight': 'normal',
         'size': 14,
         }

figure = plt.subplots(figsize=(8, 6))
plt.scatter(embedding[80:163, 1],
            embedding[80:163, 2],
            s=20,
            color='navy',
            alpha=0.5,
            label='unstress',
            )

plt.scatter(embedding[0:80, 1],
            embedding[0:80, 2],
            s=20,
            color='red',
            alpha=0.5,
            label='salt stress',
            )
plt.xlabel('Dimension 2', fontdict=font2)
plt.ylabel('Dimension 3', fontdict=font2)
plt.legend(loc="best",
           shadow=False,
           scatterpoints=1,
           prop=font3)
plt.title("Condition clustering", fontdict=font1)
plt.xticks(fontproperties='Arial',
           size=14)
plt.yticks(fontproperties='Arial',
           size=14)
plt.savefig('../singlecell/output/umap23.jpg',
            dpi=600)
plt.show()
