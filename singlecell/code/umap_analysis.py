import umap
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pandas as pd

umap_data = pd.read_excel('singlecell/output/pfba_flux.xlsx')
umap_data = umap_data.drop(0, axis=1)
umap_data.loc[:, 'metnum/20'] = umap_data.loc[:, 'metnum/20'] * 20
umap_data.loc[:, 'rxnnum/20'] = umap_data.loc[:, 'rxnnum/20'] * 20

reducer = umap.UMAP(n_components=3)
scaled_umap_data = StandardScaler().fit_transform(umap_data)
embedding = reducer.fit_transform(scaled_umap_data)


lw = 2
y = [0]*80 + [1]*83
plt.scatter(
        embedding[80:163, 0], embedding[80:163, 1], s=20, color='navy', alpha=0.3, label='unstress'
    )

plt.scatter(
        embedding[0:80, 0], embedding[0:80, 1], s=20, color='red', alpha=0.3, label='stress'
    )
plt.legend(loc="best", shadow=False, scatterpoints=1)
plt.title("umap")
#plt.xlabel('P2(2.5%)')
#plt.ylabel('P3(2.2%)')
plt.savefig('singlecell/output/umap12.jpg')
plt.show()
