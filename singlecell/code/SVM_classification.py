import pandas as pd
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import umap
from sklearn.preprocessing import StandardScaler
'''' pca
pca_data = pd.read_excel('singlecell/output/pfba_flux.xlsx')
pca_data = pca_data.drop(0, axis=1)
lable = [0]*80 + [1]*83 # 0是stress， 1是unstress

train_data, test_data = train_test_split(
    pca_data,
    random_state=1,
    train_size=0.8,
    test_size=0.2
)
train_lable, test_lable = train_test_split(
    lable,
    random_state=1,
    train_size=0.8,
    test_size=0.2
)

clf = svm.SVC(
    C=2,
    kernel='rbf',
    gamma=10,
    decision_function_shape='ovr'
)
clf.fit(train_data, train_lable)
pre_train = clf.predict(train_data)
pre_test = clf.predict(test_data)
train_acc = accuracy_score(train_lable, pre_train)
test_acc = accuracy_score(test_lable, pre_test)
all_acc = 0.8*train_acc + 0.2*test_acc
train_acc
test_acc
'''''
umap_data = pd.read_excel('singlecell/output/pfba_flux.xlsx')
umap_data = umap_data.drop(0, axis=1)
umap_data.loc[:, 'metnum/20'] = umap_data.loc[:, 'metnum/20'] * 20
umap_data.loc[:, 'rxnnum/20'] = umap_data.loc[:, 'rxnnum/20'] * 20

lable = [0]*80 + [1]*83 # 0是stress， 1是unstress
reducer = umap.UMAP(n_components=3)
scaled_umap_data = StandardScaler().fit_transform(umap_data)
embedding = reducer.fit_transform(scaled_umap_data)

train_data, test_data = train_test_split(
     umap_data,
    random_state=1,
    train_size=0.8,
    test_size=0.2
)
train_lable, test_lable = train_test_split(
    lable,
    random_state=1,
    train_size=0.8,
    test_size=0.2
)

clf = svm.SVC(
    C=2,
    kernel='rbf',
    gamma=10,
    decision_function_shape='ovr'
)
clf.fit(train_data, train_lable)
pre_train = clf.predict(train_data)
pre_test = clf.predict(test_data)
train_acc = accuracy_score(train_lable, pre_train)
test_acc = accuracy_score(test_lable, pre_test)
all_acc = 0.8*train_acc + 0.2*test_acc
all_acc