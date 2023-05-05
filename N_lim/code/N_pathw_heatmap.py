import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

# data = pd.read_excel('../output/std_pathwh.xlsx')
# data = pd.read_excel('../output/std_pathwl.xlsx')
data = pd.read_excel('../output/raw_pathw.xlsx')
idx = data.loc[:, 'idx']
data.index = data.loc[:, 'idx']
data.drop(['idx'], axis=1, inplace=True)

# sns.heatmap(data,
#             cmap="YlGnBu",
#             linewidths=1,
#             cbar_kws={"shrink": 0.8},
#             annot_kws={'family': 'Arial'})
g = sns.clustermap(data,
               figsize=(8, 5),
               # metric="correlation",
               method="complete",
               cmap="OrRd",
               linewidths=1,
               # row_cluster=False,
               #dendrogram_ratio=(.1, .2),
               cbar_pos=None
               )
# g.ax_cbar.set_title('colorbar title')
# g.ax_cbar.tick_params(axis='x', length=10)
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

plt.yticks(fontproperties='Arial',
           size=12)
plt.xticks(fontproperties='Arial',
           size=12)
plt.ylabel('')
plt.tight_layout()
# plt.savefig(
#     '../output/heatmapl.jpg',
#     dpi=600
# )
plt.savefig(
    '../output/clusterheatmap.jpg',
    dpi=600
)

plt.show()
