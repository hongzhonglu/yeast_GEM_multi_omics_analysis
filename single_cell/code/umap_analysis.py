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
    ax = plt.gca()
    ax.spines['top'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
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

reducer = umap.UMAP(n_neighbors=100, min_dist=0.5, n_components=2, metric='euclidean')
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

#umapfigure(embedding, 1, 2)



###########################################################################
# the Leiden algorithm
# import numpy as np
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import networkx as nx
# import leidenalg as la
# import igraph as ig
# from sklearn.decomposition import PCA
# from sklearn.metrics.pairwise import cosine_similarity
# import pandas as pd
# import umap
# from sklearn.metrics.pairwise import euclidean_distances
# from sklearn.manifold import TSNE
# from sklearn.preprocessing import StandardScaler
# umap_data = pd.read_excel('../output/pfba_flux.xlsx')
# umap_data = umap_data.iloc[:, 1:]
# umap_data.loc[:, 'metnum/20'] = umap_data.loc[:, 'metnum/20'] * 20
# umap_data.loc[:, 'rxnnum/20'] = umap_data.loc[:, 'rxnnum/20'] * 20
#
# similarity_matrix = cosine_similarity(umap_data)
# G = nx.Graph()
# for i in range(len(similarity_matrix)):
#     for j in range(i + 1, len(similarity_matrix)):
#         G.add_edge(i, j, weight=similarity_matrix[i, j])
#
# ig_graph = ig.Graph.from_networkx(G)
#
# partition = la.find_partition(ig_graph, la.RBConfigurationVertexPartition, resolution_parameter=1.01)
#
# labels = np.array(partition.membership)
#
# reducer = umap.UMAP(n_neighbors=100, min_dist=0.5, n_components=2, metric='euclidean')
# scaled_umap_data = StandardScaler().fit_transform(umap_data)
# reduced_data = reducer.fit_transform(scaled_umap_data)
# plt.figure(figsize=(8, 6))
# font3 = {'family': 'Arial',
#          'weight': 'normal',
#          'size': 24,
#          }
# font2 = {'family': 'Arial',
#          'weight': 'normal',
#          'size': 20,
#          }
#
# scatter = plt.scatter(reduced_data[:, 0], reduced_data[:, 1], c=labels, cmap='viridis', s=50, alpha=0.6)
# #plt.title('Leiden Clustering with UMAP Visualization', fontdict=font3)
# plt.xlabel('Component 1', fontdict=font2)
# plt.ylabel('Component 2', fontdict=font2)
# # cbar = plt.colorbar(scatter)
# # cbar.ax.tick_params(labelsize=18)
# # cbar.set_label('Cluster', fontsize=20, fontname='Arial')
# plt.xticks(fontproperties='Arial',
#            size=20)
# plt.yticks(fontproperties='Arial',
#            size=20)
# # for label in cbar.ax.get_yticklabels():
# #     label.set_fontname('Arial')
#
# plt.tight_layout()
# plt.savefig('./Leiden.tif')
# plt.show()
#
# ##############################
# import numpy as np
# import matplotlib.pyplot as plt
# from sklearn.cluster import AgglomerativeClustering
# from sklearn.manifold import TSNE
# from sklearn.preprocessing import StandardScaler
#
# scaler = StandardScaler()
# data_scaled = scaler.fit_transform(umap_data)
#
# n_clusters = 2
# clustering = AgglomerativeClustering(n_clusters=n_clusters, metric='euclidean', linkage='ward')
# labels = clustering.fit_predict(data_scaled)
#
# print("社区标签分布:", np.unique(labels, return_counts=True))
# reducer = umap.UMAP(n_neighbors=100, min_dist=0.5, n_components=2, metric='euclidean')
# scaled_umap_data = StandardScaler().fit_transform(umap_data)
# reduced_data = reducer.fit_transform(scaled_umap_data)
# plt.figure(figsize=(8, 6))
# font3 = {'family': 'Arial',
#          'weight': 'normal',
#          'size': 24,
#          }
# font2 = {'family': 'Arial',
#          'weight': 'normal',
#          'size': 20,
#          }
#
# scatter = plt.scatter(reduced_data[:, 0], reduced_data[:, 1], c=labels, cmap='viridis', s=50, alpha=0.6)
# #plt.title('Hierarchical Clustering with UMAP Visualization', fontdict=font3)
# plt.xlabel('Component 1', fontdict=font2)
# plt.ylabel('Component 2', fontdict=font2)
# #cbar = plt.colorbar(scatter)
# # cbar.ax.tick_params(labelsize=18)
# # cbar.set_label('Cluster', fontsize=20, fontname='Arial')
# plt.xticks(fontproperties='Arial',
#            size=20)
# plt.yticks(fontproperties='Arial',
#            size=20)
# # for label in cbar.ax.get_yticklabels():
# #     label.set_fontname('Arial')
#
# plt.tight_layout()
# plt.savefig('./Hierarchical.tif')
# plt.show()