##################################################
import pandas as pd
from sklearn.metrics import accuracy_score
import umap
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.decomposition import KernelPCA


def rf_kpca(data, lable, method, datatype):
    train_data, test_data, train_lable, test_lable = train_test_split(
        data,
        lable,
        random_state=1,
        train_size=0.7,
        test_size=0.3,
    )

    pipeline = Pipeline([
        ('kpca', KernelPCA()),
        ('rf', RandomForestClassifier(random_state=42))
    ])

    # 定义要搜索的参数网格
    param_grid = {
        'kpca__n_components': [5, 10, 100, 500],
        'kpca__kernel': ['linear', 'poly', 'rbf'],
        'rf__n_estimators': [10, 100, 500],
        'rf__max_depth': [None, 10, 20],
        'rf__min_samples_split': [2, 5, 10],
        'rf__min_samples_leaf': [1, 2, 4],
    }

    #grid_search = GridSearchCV(pipeline, param_grid, cv=5)
    grid_search = GridSearchCV(pipeline, param_grid, cv=5, n_jobs=-1)

    grid_search.fit(train_data, train_lable)
    best_model = grid_search.best_estimator_
    pre_train = best_model.predict(train_data)
    pre_test = best_model.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)
    print("Flux test_acc：", test_acc)
    print("Flux train_acc：", train_acc)
    ss = 0  # exp stress-simu stress
    s_ns = 0  # stress-unstress
    ns_s = 0
    ns_ns = 0
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

    print("Best parameters found: ", grid_search.best_params_)
    print("Best cross-validation score: {:.2f}".format(grid_search.best_score_))
    result = pd.DataFrame(index=['value_{}_{}'.format(method, datatype)], columns=['test_acc', 'train_acc', 'ss', 's_ns', 'ns_ns', 'ns_s'])
    result.loc['value_{}_{}'.format(method, datatype), 'test_acc'] = test_acc
    result.loc['value_{}_{}'.format(method, datatype), 'train_acc'] = train_acc
    result.loc['value_{}_{}'.format(method, datatype), 'ss'] = ss
    result.loc['value_{}_{}'.format(method, datatype), 's_ns'] = s_ns
    result.loc['value_{}_{}'.format(method, datatype), 'ns_ns'] = ns_ns
    result.loc['value_{}_{}'.format(method, datatype), 'ns_s'] = ns_s

    return result


def rf(data, lable, method, datatype):
    train_data, test_data, train_lable, test_lable = train_test_split(
        data,
        lable,
        random_state=1,
        train_size=0.7,
        test_size=0.3,
    )
    param_grid = {
        'n_estimators': [10, 100, 500],
        'max_depth': [None, 10, 20],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4],
    }
    rf_classifier = RandomForestClassifier(random_state=42)

    grid_search = GridSearchCV(rf_classifier, param_grid, cv=5, n_jobs=-1)
    #grid_search = GridSearchCV(rf_classifier, param_grid, cv=5)

    grid_search.fit(train_data, train_lable)

    print("Best Parameters:", grid_search.best_params_)

    best_rf_classifier = RandomForestClassifier(**grid_search.best_params_, random_state=42)

    best_rf_classifier.fit(train_data, train_lable)

    y_pred = best_rf_classifier.predict(test_data)

    pre_train = best_rf_classifier.predict(train_data)
    pre_test = best_rf_classifier.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)

    ss = 0  # exp stress-simu stress
    s_ns = 0  # stress-unstress
    ns_s = 0
    ns_ns = 0
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

    result = pd.DataFrame(index=['value_{}_{}'.format(method, datatype)], columns=['test_acc', 'train_acc', 'ss', 's_ns', 'ns_ns', 'ns_s'])
    result.loc['value_{}_{}'.format(method, datatype), 'test_acc'] = test_acc
    result.loc['value_{}_{}'.format(method, datatype), 'train_acc'] = train_acc
    result.loc['value_{}_{}'.format(method, datatype), 'ss'] = ss
    result.loc['value_{}_{}'.format(method, datatype), 's_ns'] = s_ns
    result.loc['value_{}_{}'.format(method, datatype), 'ns_ns'] = ns_ns
    result.loc['value_{}_{}'.format(method, datatype), 'ns_s'] = ns_s
    return result


data = pd.read_excel('../output/pfba_flux.xlsx')
data = data.iloc[:, 1:]
data.loc[:, 'metnum/20'] = data.loc[:, 'metnum/20'] * 20
data.loc[:, 'rxnnum/20'] = data.loc[:, 'rxnnum/20'] * 20

lable = [0] * 80 + [1] * 83  # 0 is stress， 1 is unstress
reducer = umap.UMAP(n_components=5)
scaled_umap_flux = StandardScaler().fit_transform(data)
embedding_flux = reducer.fit_transform(scaled_umap_flux)
transdata = pd.read_excel('../data/GSE102475_GASCH_NaCl-scRNAseq_NormData.xlsx', index_col=0)
transdata = transdata.T


print('rf_umap flux:')
rf_umap_flux = rf(embedding_flux, lable, 'umap', 'flux')
rf_umap_flux.to_excel('../output/rf_umap_flux.xlsx')

print('rf_umap trans:')
scaled_umap_trans = StandardScaler().fit_transform(transdata)
embedding_trans = reducer.fit_transform(scaled_umap_trans)
rf_umap_trans = rf(embedding_trans, lable, 'umap', 'trans')
rf_umap_trans.to_excel('../output/rf_umap_trans.xlsx')

print('rf flux:')
rf_flux = rf(data, lable, 'off', 'flux')
rf_flux.to_excel('../output/rf_flux.xlsx')

print('rf trans:')
rf_trans = rf(transdata, lable, 'off', 'trans')
rf_trans.to_excel('../output/rf_trans.xlsx')

print('rf_kpca flux:')
rf_kpca_flux = rf_kpca(data, lable, 'kpca', 'flux')
rf_kpca_flux.to_excel('../output/rf_kpca_flux.xlsx')
print('rf_kpca trans')
rf_kpca_trans = rf_kpca(transdata, lable, 'kpca', 'trans')
rf_kpca_trans.to_excel('../output/rf_kpca_trans.xlsx')


result = pd.concat([rf_kpca_flux, rf_kpca_trans, rf_umap_flux,
                    rf_umap_trans, rf_flux, rf_trans], axis=0)

result.to_excel('../output/result.xlsx')

#############################################################
# randomization of cell labels
# import random
# lable = [0] * 80 + [1] * 83  # 0 is stress， 1 is unstress
# for i in range(10):
#     shuffled_list = random.sample(lable, len(lable))
#     print(i)
#     rf_kpca_flux = rf_kpca(data, shuffled_list, 'kpca', 'flux')


