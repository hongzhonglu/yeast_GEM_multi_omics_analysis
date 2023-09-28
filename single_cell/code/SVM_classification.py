##################################################
# apply SVM to classify two condition
from sklearn.naive_bayes import GaussianNB
import pandas as pd
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import umap
from sklearn.preprocessing import StandardScaler
import os

os.chdir('..')

umap_data = pd.read_excel('../singlecell/output/pfba_flux.xlsx')
umap_data = umap_data.drop(0, axis=1)
umap_data.loc[:, 'metnum/20'] = umap_data.loc[:, 'metnum/20'] * 20
umap_data.loc[:, 'rxnnum/20'] = umap_data.loc[:, 'rxnnum/20'] * 20

lable = [0] * 80 + [1] * 83  # 0 is stress， 1 is unstress
reducer = umap.UMAP(n_components=3)
scaled_umap_data = StandardScaler().fit_transform(umap_data)
embedding = reducer.fit_transform(scaled_umap_data)
# train data size : test data size = 8:2
# train_data, test_data = train_test_split(
#     umap_data,
#     random_state=1,
#     train_size=0.8,
#     test_size=0.2
# )
# train_lable, test_lable = train_test_split(
#     lable,
#     random_state=1,
#     train_size=0.8,
#     test_size=0.2
# )
#
# clf = svm.SVC(
#     C=2,
#     kernel='rbf',
#     gamma=10,
#     decision_function_shape='ovr'
# )
# clf.fit(train_data, train_lable)
# pre_train = clf.predict(train_data)
# pre_test = clf.predict(test_data)
# train_acc = accuracy_score(train_lable, pre_train)
# test_acc = accuracy_score(test_lable, pre_test)
train_data, test_data, train_lable, test_lable = train_test_split(
    umap_data,
    lable,
    random_state=1,
    train_size=0.8,
    test_size=0.2,
)
gnb = GaussianNB()
gnb.fit(train_data, train_lable)

pre_train = gnb.predict(train_data)
pre_test = gnb.predict(test_data)
train_acc = accuracy_score(train_lable, pre_train)
test_acc = accuracy_score(test_lable, pre_test)


ss = 0 #exp stress-simu stress
s_ns = 0 #stress-unstress
ns_s = 0
ns_ns = 0
fp_l = []
fn_l = []
for ff in range(len(pre_train)):
    if train_lable[ff] == 0 and pre_train[ff] == 0:
        ss += 1
    if train_lable[ff] == 1 and pre_train[ff] == 1:
        ns_ns += 1
    if train_lable[ff] == 0 and pre_train[ff] == 1:
        s_ns += 1
    if train_lable[ff] == 1 and pre_train[ff] == 0:
        ns_s += 1
for ff in range(len(pre_test)):
    if test_lable[ff] == 0 and pre_test[ff] == 0:
        ss += 1
    if test_lable[ff] == 1 and pre_test[ff] == 1:
        ns_ns += 1
    if test_lable[ff] == 0 and pre_test[ff] == 1:
        s_ns += 1
    if test_lable[ff] == 1 and pre_test[ff] == 0:
        ns_s += 1

print(
    'accuracy:' + str(test_acc) + '\n' +
    'ss:' + str(ss) + '\n' +
    's_ns:' + str(s_ns) + '\n' +
    'ns_ns:' + str(ns_ns) + '\n' +
    'ns_s:' + str(ns_s) + '\n'
)
