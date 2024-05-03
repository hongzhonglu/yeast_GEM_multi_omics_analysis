############################################
# apply umap to visualize classification
# get "umap.jpg"
# first run "get_pfba_flux.py" to compute flux data "pfba_flux.xlsx"

import umap
from sklearn.preprocessing import StandardScaler
import pandas as pd
import matplotlib.pyplot as plt

def umapfigure(embedding, a, b):
    font3 = {'family': 'Arial',
             'weight': 'normal',
             'size': 26,
             }

    plt.subplots(figsize=(8, 6))
    plt.scatter(embedding[80:163, a - 1],
                embedding[80:163, b - 1],
                s=20,
                color='#3174A1',
                alpha=0.6,
                label='Unstressed',
                )

    plt.scatter(embedding[0:80, a - 1],
                embedding[0:80, b - 1],
                s=20,
                color='#E1812C',
                alpha=0.6,
                label='Stress',
                )
    plt.xlabel('Dimension {}'.format(str(a)), fontdict=font3)
    plt.ylabel('Dimension {}'.format(str(b)), fontdict=font3)
    plt.legend(loc="best",
               shadow=False,
               scatterpoints=1,
               prop=font3)
    plt.xticks(fontproperties='Arial',
               size=24)
    plt.yticks(fontproperties='Arial',
               size=26)
    plt.tight_layout()
    plt.savefig('../output/umap{}{}.jpg'.format(str(a), str(b)),
                dpi=600)
    plt.show()


umap_data = pd.read_excel('../output/pfba_flux.xlsx')
umap_data = umap_data.iloc[:, 1:]
umap_data.loc[:, 'metnum/20'] = umap_data.loc[:, 'metnum/20'] * 20
umap_data.loc[:, 'rxnnum/20'] = umap_data.loc[:, 'rxnnum/20'] * 20

reducer = umap.UMAP(n_neighbors=100, min_dist=0.8, n_components=2, metric='euclidean')
scaled_umap_data = StandardScaler().fit_transform(umap_data)
embedding = reducer.fit_transform(scaled_umap_data)

# Testing the independence of features
# df = pd.DataFrame(embedding)
# features = list(df.columns)
# for i in range(len(features)):
#     for j in range(i + 1, len(features)):
#         contingency_table = pd.crosstab(df[features[i]], df[features[j]])
#         chi2, p, dof, expected = chi2_contingency(contingency_table)
#
#         print(f"result - {features[i]} vs {features[j]}:")
#         print(f"K value: {chi2}, p: {p}\n")

umapfigure(embedding, 1, 2)
