
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

pca_data = pd.read_excel('singlecell/output/pfba_flux.xlsx')
pca_data = pca_data.drop(0, axis=1)
pca_data = pca_data.drop('metnum/20', axis=1)
pca_data = pca_data.drop('rxnnum/20', axis=1)
pca = PCA(n_components=5)
X_r = pca.fit(pca_data).transform(pca_data)
Y_r = pca.fit(pca_data.T).transform(pca_data.T)
Y_df = pd.DataFrame(Y_r, index=pca_data.columns.values.tolist())
