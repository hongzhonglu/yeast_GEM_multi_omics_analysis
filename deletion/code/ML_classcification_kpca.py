import pandas as pd
from sklearn import svm
from sklearn.metrics import accuracy_score, recall_score, classification_report
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
import logging
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.decomposition import KernelPCA
from datetime import datetime
from tqdm import tqdm


logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')


def kpca_reduce(ml_data):
    # components = [50, 100, 500, 1000]
    # kernel = ['linear', 'poly', 'rbf']
    #########################################
    # trans
    # gnb/knn best
    # components = [1000]
    # kernel = ['rbf']
    # # log/mlp
    # components = [500]
    # kernel = ['linear']
    # rf/svm
    # components = [50]
    # kernel = ['linear']
    #########################################
    # flux
    # gnb best
    components = [500]
    kernel = ['rbf']
    # knn best
    # components = [100]
    # kernel = ['rbf']
    # log/mlp/svm
    # components = [500]
    # kernel = ['linear']
    # rf
    # components = [50]
    # kernel = ['linear']
    #########################################
    # ft
    # gnb/rf best
    # components = [500]
    # kernel = ['rbf']
    # knn/log best
    # components = [50]
    # kernel = ['rbf']
    # mlp
    # components = [500]
    # kernel = ['poly']
    # svm
    # components = [100]
    # kernel = ['rbf']

    X_transformed = {}
    for c in components:
        for k in kernel:
            transformer = KernelPCA(n_components=c, kernel=k)
            X_transformed[str(c) + '' + k] = transformer.fit_transform(ml_data)
    return X_transformed


def svm_class(X_transformed, lable, data):
    final = pd.DataFrame()
    n = 0
    for k, ml in X_transformed.items():
        train_data, test_data, train_lable, test_lable = train_test_split(
            ml,
            lable,
            random_state=1,
            train_size=0.7,
            test_size=0.3,
        )
        clf = svm.SVC(probability=True)
        param_grid = {
            'C': [0.1, 1, 10, 100],
            'gamma': [0.1, 0.01, 0.001],
            # 'linear' and 'poly' cause that the simulation time is too long.
            'kernel': ['rbf', 'sigmoid']
        }
        # t best
        # param_grid = {
        #     'C': [100],
        #     'gamma': [0.001],
        #     # 'linear' and 'poly' cause that the simulation time is too long.
        #     'kernel': ['rbf']
        # }
        # f best
        # param_grid = {
        #     'C': [100],
        #     'gamma': [0.1],
        #     # 'linear' and 'poly' cause that the simulation time is too long.
        #     'kernel': ['rbf']
        # }
        # ft best
        # param_grid = {
        #     'C': [100],
        #     'gamma': [0.1],
        #     # 'linear' and 'poly' cause that the simulation time is too long.
        #     'kernel': ['rbf']
        # }
        grid_search = GridSearchCV(clf, param_grid, cv=5, verbose=3, scoring='accuracy', n_jobs=-1)
        grid_search.fit(train_data, train_lable)

        pre_train_pro = grid_search.predict_proba(train_data)
        pre_test_pro = grid_search.predict_proba(test_data)
        pre_pro_svm_ft = np.vstack((pre_train_pro, pre_test_pro))
        pre_train = grid_search.predict(train_data)
        pre_test = grid_search.predict(test_data)
        true_data_svm_ft = np.hstack((train_lable, test_lable))


        final.loc[n, 'train_acc'] = accuracy_score(train_lable, pre_train)
        final.loc[n, 'test_acc'] = accuracy_score(test_lable, pre_test)
        final.loc[n, 'recall_t'] = recall_score(test_lable, pre_test, average="macro")
        final.loc[n, 'fscore_t'] = f1_score(test_lable, pre_test, average="macro")
        final.loc[n, 'best_params'] = str(grid_search.best_params_)
        final.loc[n, 'kpca'] = k

        n += 1
    final.sort_values(by='test_acc', inplace=True, ascending=False)
    best = final.iloc[0, :]
    current_time = datetime.now().strftime("%H h %M min")
    best.to_excel('../output/kpca_Reduction/svm_best_{}.xlsx'.format(data))
    return best


def MLP_class(X_transformed, lable, data):
    final = pd.DataFrame()
    n = 0
    for k, ml in X_transformed.items():
        train_data, test_data, train_lable, test_lable = train_test_split(
            ml,
            lable,
            random_state=1,
            train_size=0.7,
            test_size=0.3,
        )
        clf = MLPClassifier()
        param_grid = {
            'hidden_layer_sizes': [(10,), (10, 10), [(10, 20, 10)]],
            'alpha': [0.0001, 0.001, 0.01],
            'activation': ['logistic', 'tanh', 'relu']
        }
        # t
        # param_grid = {
        #     'hidden_layer_sizes': [(10,)],
        #     'alpha': [0.01],
        #     'activation': ['tanh']
        # }
        # f
        # param_grid = {
        #     'hidden_layer_sizes': [(10, 10)],
        #     'alpha': [0.0001],
        #     'activation': ['relu']
        # }
        # ft
        # param_grid = {
        #     'hidden_layer_sizes': [(10, 10)],
        #     'alpha': [0.0001],
        #     'activation': ['logistic']
        # }

        grid_search = GridSearchCV(clf, param_grid, cv=5, verbose=2, scoring='accuracy')
        grid_search.fit(train_data, train_lable)

        pre_train_pro = grid_search.predict_proba(train_data)
        pre_test_pro = grid_search.predict_proba(test_data)
        pre_pro_mlp_ft = np.vstack((pre_train_pro, pre_test_pro))
        pre_train = grid_search.predict(train_data)
        pre_test = grid_search.predict(test_data)
        true_data_mlp_ft = np.hstack((train_lable, test_lable))

        final.loc[n, 'train_acc'] = accuracy_score(train_lable, pre_train)
        final.loc[n, 'test_acc'] = accuracy_score(test_lable, pre_test)
        final.loc[n, 'recall_t'] = recall_score(test_lable, pre_test, average="macro")
        final.loc[n, 'fscore_t'] = f1_score(test_lable, pre_test, average="macro")
        final.loc[n, 'best_params'] = str(grid_search.best_params_)
        final.loc[n, 'kpca'] = k

        n += 1
    final.sort_values(by='test_acc', inplace=True, ascending=False)
    best = final.iloc[0, :]
    current_time = datetime.now().strftime("%H h %M min")
    best.to_excel('../output/kpca_Reduction/mlp_best_{}.xlsx'.format(data))
    return best


def bayes_class(X_transformed, lable, data):
    final = pd.DataFrame()
    n = 0
    for k, ml in X_transformed.items():
        train_data, test_data, train_lable, test_lable = train_test_split(
            ml,
            lable,
            random_state=1,
            train_size=0.7,
            test_size=0.3,
        )
        gnb = GaussianNB()
        # param_grid = {'var_smoothing': [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]}
        # t
        #param_grid = {'var_smoothing': [1e-9]}
        # f/ft
        param_grid = {'var_smoothing': [1e-5]}


        grid_search = GridSearchCV(gnb, param_grid, cv=2, verbose=2, scoring='accuracy', n_jobs=-1)
        grid_search.fit(train_data, train_lable)

        pre_train_pro = grid_search.predict_proba(train_data)
        pre_test_pro = grid_search.predict_proba(test_data)
        pre_pro_gnb_ft = np.vstack((pre_train_pro, pre_test_pro))
        pre_train = grid_search.predict(train_data)
        pre_test = grid_search.predict(test_data)
        true_data_gnb_ft = np.hstack((train_lable, test_lable))


        final.loc[n, 'train_acc'] = accuracy_score(train_lable, pre_train)
        final.loc[n, 'test_acc'] = accuracy_score(test_lable, pre_test)
        final.loc[n, 'recall_t'] = recall_score(test_lable, pre_test, average="macro")
        final.loc[n, 'fscore_t'] = f1_score(test_lable, pre_test, average="macro")
        final.loc[n, 'best_params'] = str(grid_search.best_params_)
        final.loc[n, 'kpca'] = k
        n += 1
    final.sort_values(by='test_acc', inplace=True, ascending=False)
    best = final.iloc[0, :]
    current_time = datetime.now().strftime("%H h %M min")
    best.to_excel('../output/kpca_Reduction/gnb_best_{}.xlsx'.format(data))
    return best


def randomf_class(X_transformed, lable, data):
    final = pd.DataFrame()
    n = 0
    for k, ml in X_transformed.items():
        train_data, test_data, train_lable, test_lable = train_test_split(
            ml,
            lable,
            random_state=1,
            train_size=0.7,
            test_size=0.3,
        )
        rf = RandomForestClassifier()
        param_grid = {
            'n_estimators': [10, 100, 300],
            'max_depth': [None, 5, 10],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 4],
        }
        # t
        # param_grid = {
        #     'n_estimators': [100],
        #     'max_depth': [10],
        #     'min_samples_split': [2],
        #     'min_samples_leaf': [1],
        # }
        # f
        # param_grid = {
        #     'n_estimators': [300],
        #     'max_depth': [None],
        #     'min_samples_split': [2],
        #     'min_samples_leaf': [2],
        # }
        # ft
        # param_grid = {
        #     'n_estimators': [300],
        #     'max_depth': [None],
        #     'min_samples_split': [5],
        #     'min_samples_leaf': [1],
        # }
        grid_search = GridSearchCV(rf, param_grid, cv=5, verbose=2, scoring='accuracy', n_jobs=-1)
        grid_search.fit(train_data, train_lable)

        pre_train_pro = grid_search.predict_proba(train_data)
        pre_test_pro = grid_search.predict_proba(test_data)
        pre_pro = np.vstack((pre_train_pro, pre_test_pro))
        pre_pro_rf_ft = np.vstack((pre_train_pro, pre_test_pro))
        pre_train = grid_search.predict(train_data)
        pre_test = grid_search.predict(test_data)
        true_data = np.hstack((train_lable, test_lable))
        true_data_rf_ft = np.hstack((train_lable, test_lable))

        final.loc[n, 'train_acc'] = accuracy_score(train_lable, pre_train)
        final.loc[n, 'test_acc'] = accuracy_score(test_lable, pre_test)
        final.loc[n, 'recall_t'] = recall_score(test_lable, pre_test, average="macro")
        final.loc[n, 'fscore_t'] = f1_score(test_lable, pre_test, average="macro")
        final.loc[n, 'best_params'] = str(grid_search.best_params_)
        final.loc[n, 'kpca'] = k
        n += 1
    final.sort_values(by='test_acc', inplace=True, ascending=False)
    best = final.iloc[0, :]
    current_time = datetime.now().strftime("%H h %M min")
    best.to_excel('../output/kpca_Reduction/rf_best_{}.xlsx'.format(data))
    return best


def knn_class(X_transformed, lable, data):
    final = pd.DataFrame()
    n = 0
    for k, ml in X_transformed.items():
        train_data, test_data, train_lable, test_lable = train_test_split(
            ml,
            lable,
            random_state=1,
            train_size=0.7,
            test_size=0.3,
        )
        param_grid = {
            'n_neighbors': [3, 5, 7, 9],
            'weights': ['uniform', 'distance'],
            'metric': ['euclidean', 'manhattan', 'minkowski'],
            'algorithm': ['auto', 'ball_tree', 'kd_tree']
        }
        # t
        # param_grid = {
        #     'n_neighbors': [5],
        #     'weights': ['distance'],
        #     'metric': ['manhattan'],
        #     'algorithm': ['auto']
        # }
        # f
        # param_grid = {
        #     'n_neighbors': [3],
        #     'weights': ['distance'],
        #     'metric': ['euclidean'],
        #     'algorithm': ['auto']
        # }
        # ft
        # param_grid = {
        #     'n_neighbors': [3],
        #     'weights': ['distance'],
        #     'metric': ['euclidean'],
        #     'algorithm': ['auto']
        # }
        knn = KNeighborsClassifier()
        grid_search = GridSearchCV(knn, param_grid, cv=5, scoring='accuracy', verbose=2, n_jobs=-1)
        grid_search.fit(train_data, train_lable)

        pre_train_pro = grid_search.predict_proba(train_data)
        pre_test_pro = grid_search.predict_proba(test_data)
        pre_pro = np.vstack((pre_train_pro, pre_test_pro))
        pre_pro_knn_ft = np.vstack((pre_train_pro, pre_test_pro))
        pre_train = grid_search.predict(train_data)
        pre_test = grid_search.predict(test_data)
        true_data = np.hstack((train_lable, test_lable))
        true_data_knn_ft = np.hstack((train_lable, test_lable))

        final.loc[n, 'train_acc'] = accuracy_score(train_lable, pre_train)
        final.loc[n, 'test_acc'] = accuracy_score(test_lable, pre_test)
        final.loc[n, 'recall_t'] = recall_score(test_lable, pre_test, average="macro")
        final.loc[n, 'fscore_t'] = f1_score(test_lable, pre_test, average="macro")
        final.loc[n, 'best_params'] = str(grid_search.best_params_)
        final.loc[n, 'kpca'] = k

        n += 1
    final.sort_values(by='test_acc', inplace=True, ascending=False)
    best = final.iloc[0, :]
    current_time = datetime.now().strftime("%H h %M min")
    best.to_excel('../output/kpca_Reduction/knn_best_{}.xlsx'.format(data))
    return best


def logistic_class(X_transformed, lable, data):
    final = pd.DataFrame()
    n = 0
    for k, ml in X_transformed.items():
        train_data, test_data, train_lable, test_lable = train_test_split(
            ml,
            lable,
            random_state=1,
            train_size=0.7,
            test_size=0.3,
        )
        param_grid = {
            'C': [0.01, 0.1, 1],
            'solver': ['liblinear', 'saga'],
            'class_weight': [None, 'balanced']
        }
        # t
        # param_grid = {
        #     'C': [1],
        #     'solver': ['liblinear'],
        #     'class_weight': [None]
        # }
        # f
        # param_grid = {
        #     'C': [1],
        #     'solver': ['liblinear'],
        #     'class_weight': ['balanced']
        # }
        # ft
        # param_grid = {
        #     'C': [1],
        #     'solver': ['saga'],
        #     'class_weight': ['balanced']
        # }

        clf = LogisticRegression()
        grid_search = GridSearchCV(clf, param_grid, cv=5, scoring='accuracy', verbose=2, n_jobs=-1)
        grid_search.fit(train_data, train_lable)

        pre_train_pro = grid_search.predict_proba(train_data)
        pre_test_pro = grid_search.predict_proba(test_data)
        pre_pro_log_ft = np.vstack((pre_train_pro, pre_test_pro))
        pre_train = grid_search.predict(train_data)
        pre_test = grid_search.predict(test_data)
        true_data_log_ft = np.hstack((train_lable, test_lable))

        final.loc[n, 'train_acc'] = accuracy_score(train_lable, pre_train)
        final.loc[n, 'test_acc'] = accuracy_score(test_lable, pre_test)
        final.loc[n, 'recall_t'] = recall_score(test_lable, pre_test, average="macro")
        final.loc[n, 'fscore_t'] = f1_score(test_lable, pre_test, average="macro")
        final.loc[n, 'best_params'] = str(grid_search.best_params_)
        final.loc[n, 'kpca'] = k

        n += 1
    final.sort_values(by='test_acc', inplace=True, ascending=False)
    best = final.iloc[0, :]
    current_time = datetime.now().strftime("%H h %M min")
    best.to_excel('../output/kpca_Reduction/log_best_{}.xlsx'.format(data))
    return best


def prepare_mldata(sortedgene, flux_or_trans):
    ml_data = pd.DataFrame(index=sortedgene.index, columns=flux_or_trans.columns)
    for g in sortedgene.index:
        ml_data.loc[g, :] = flux_or_trans.loc[g, :]
    ml_data.index = range(len(ml_data.index))
    return ml_data


def roc(true_data, pre_pro, n, method, reduction, datatype):
    if n == 50:
        n_class = 9
    elif n == 30:
        n_class = 13
    elif n == 90:
        n_class = 4
    elif n == 70:
        n_class = 5
    lable = label_binarize(true_data, np.arange(n_class))  # 装换成类似二进制的编码
    fpr, tpr, thresholds = metrics.roc_curve(lable.ravel(), pre_pro.ravel())
    auc = metrics.auc(fpr, tpr)
    # 绘图
    font2 = {'family': 'Arial',
             'weight': 'normal',
             'size': 14,
             }
    font3 = {'family': 'Arial',
             'weight': 'normal',
             'size': 12,
             }
    mpl.rcParams['font.sans-serif'] = u'SimHei'
    mpl.rcParams['axes.unicode_minus'] = False
    # FPR就是横坐标,TPR就是纵坐标
    plt.plot(fpr, tpr, c='r', lw=2, alpha=0.7, label=u'AUC=%.3f' % auc)
    plt.plot((0, 1), (0, 1), c='#808080', lw=1, ls='--', alpha=0.7)
    plt.xlim((-0.01, 1.02))
    plt.ylim((-0.01, 1.02))
    plt.xticks(np.arange(0, 1.1, 0.1), fontproperties='Arial', size=12)
    plt.yticks(np.arange(0, 1.1, 0.1), fontproperties='Arial', size=12)
    plt.xlabel('False Positive Rate', fontdict=font2)
    plt.ylabel('True Positive Rate', fontdict=font2)
    plt.grid(b=True, ls=':')
    plt.legend(loc='lower right', fancybox=True, framealpha=0.8, prop=font3)
    plt.title('ROC and AUC of {}'.format(method), fontdict=font2)
    plt.tight_layout()
    plt.savefig('../output/kpca_Reduction/ROC and AUC of {} {} {} {}.png'.format(method, n, reduction, datatype),
                dpi=600,
                transparent=True)
    plt.show()
    logging.info('figure saved')


def run_ml(ml_data, lable, n, data):
    X_transformed = kpca_reduce(ml_data)
    try:
        svmbest = svm_class(X_transformed, lable, data)
    except:
        print('svm error')
    try:
        MLPbest = MLP_class(X_transformed, lable, data)
    except:
        print('MLP error')
    try:
        bayesbest = bayes_class(X_transformed, lable, data)
    except:
        print('bayes error')
    try:
        ranfbest = randomf_class(X_transformed, lable, data)
    except:
        print('bayes error')
    try:
        knnbest = knn_class(X_transformed, lable, data)
    except:
        print('bayes error')
    try:
        logisbest = logistic_class(X_transformed, lable, data)
    except:
        print('bayes error')
    # output = pd.DataFrame(index=['svm_{}'.format(n), 'MLP_{}'.format(n), 'bayes_{}'.format(n),
    #                              'randomforest_{}'.format(n), 'knn_{}'.format(n), 'logistic_{}'.format(n)],
    #                       columns=['true_data', 'pre_pro'])
    # output.loc['svm_{}'.format(n), 'true_data'] = [svmtrue_data, 'kpca', svmprepro]
    # output.loc['MLP_{}'.format(n), 'true_data'] = [MLPtrue_data, 'kpca', MLPprepro]
    # output.loc['bayes_{}'.format(n), 'true_data'] = [bayestrue_data, 'kpca', bayesprepro]
    # output.loc['randomforest_{}'.format(n), 'true_data'] = [ranftrue_data, 'kpca', ranfprepro]
    # output.loc['knn_{}'.format(n), 'true_data'] = [knntrue_data, 'kpca', knnprepro]
    # output.loc['logistic_{}'.format(n), 'true_data'] = [logistrue_data, 'kpca', logisprepro]
    #
    # output.loc['svm_{}'.format(n), 'pre_pro'] = [svmtrue_data, 'kpca', svmprepro]
    # output.loc['MLP_{}'.format(n), 'pre_pro'] = [MLPtrue_data, 'kpca', MLPprepro]
    # output.loc['bayes_{}'.format(n), 'pre_pro'] = [bayestrue_data, 'kpca', bayesprepro]
    # output.loc['randomforest_{}'.format(n), 'pre_pro'] = [ranftrue_data, 'kpca', ranfprepro]
    # output.loc['knn_{}'.format(n), :] = [knntrue_data, 'kpca', knnprepro]
    # output.loc['logistic_{}'.format(n), 'pre_pro'] = [logistrue_data, 'kpca', logisprepro]

    #return output


def prepare_gene_class(n):
    geneclass = pd.read_excel('../data/geneclass.xlsx')
    geneclass.index = geneclass.loc[:, 'gene']
    if n == 50:
        sortedgene = pd.read_excel('../data/SortedFluxGene_50.xlsx')
        sortedgene.index = sortedgene.loc[:, 'genes']
        for s in sortedgene.index:
            sortedgene.loc[s, 'orf name'] = geneclass.loc[s, 'orf name']
        lable = [0] * 69 + [1] * 57 + [3] * 72 + [4] * 137 \
                + [5] * 148 + [6] * 92 + [7] * 64 + [8] * 50 \
                + [9] * 110
    elif n == 90:
        sortedgene = pd.read_excel('../data/SortedFluxGene_90.xlsx')
        sortedgene.index = sortedgene.loc[:, 'genes']
        for s in sortedgene.index:
            sortedgene.loc[s, 'orf name'] = geneclass.loc[s, 'orf name']
        lable = [0] * 137 + [1] * 148 + [3] * 92 + [4] * 110
    elif n == 30:
        sortedgene = pd.read_excel('../data/SortedFluxGene_30.xlsx')
        sortedgene.index = sortedgene.loc[:, 'genes']
        for s in sortedgene.index:
            sortedgene.loc[s, 'orf name'] = geneclass.loc[s, 'orf name']
        lable = [0] * 69 + [1] * 57 + [3] * 72 + [4] * 137 \
                + [5] * 148 + [6] * 92 + [7] * 64 + [8] * 50 \
                + [9] * 110 + [10] * 44 + [11] * 35 + [12] * 34 \
                + [13] * 35
    elif n == 70:
        sortedgene = pd.read_excel('../data/SortedFluxGene_70.xlsx')
        sortedgene.index = sortedgene.loc[:, 'genes']
        for s in sortedgene.index:
            sortedgene.loc[s, 'orf name'] = geneclass.loc[s, 'orf name']
        lable = [0] * 72 + [1] * 137 \
                + [2] * 148 + [3] * 92 \
                + [4] * 110
    return sortedgene, lable


def get_result(n):
    # prepare flux, boounds and trans
    msb = pd.read_csv('../data/msbdata.csv')
    flux = pd.read_csv('../output/Glucose_fluxDataset.csv')
    sortedgene, lable = prepare_gene_class(n)
    trans_data = msb.drop(columns=['geneName', 'commonName']).T
    flux_data = flux.drop(columns=['Var1']).T
    ml_trans_data = prepare_mldata(sortedgene, trans_data)
    ml_flux_data = prepare_mldata(sortedgene, flux_data)
    ft = pd.concat([ml_trans_data, ml_flux_data], axis=1)
    logging.info('data loaded')

    # only with trans
    only_trans = run_ml(ml_trans_data, lable, n)
    # for m in only_trans.index.tolist():
    #     roc(only_trans.loc[m, 'true_data'], only_trans.loc[m, 'pre_pro'], n, m, only_trans.loc[m, 'reduced'],
    #         'only_trans')
    logging.info('trans finished')
    # only with flux
    only_flux = run_ml(ml_flux_data, lable, n)
    # for m in only_flux.index.tolist():
    #     roc(only_flux.loc[m, 'true_data'], only_flux.loc[m, 'pre_pro'], n, m, only_flux.loc[m, 'reduced'], 'only_flux')
    # logging.info('flux finished')
    # flux and trans
    trans_flux = run_ml(ft, lable, n)
    # for m in only_flux.index.tolist():
    #     roc(trans_flux.loc[m, 'true_data'], trans_flux.loc[m, 'pre_pro'], n, m, trans_flux.loc[m, 'reduced'],
    #         'trans_flux')
    logging.info('all finished')
    return only_trans, only_flux, trans_flux


def from_result_to_output():
    only_trans_data = pd.DataFrame()
    only_flux_data = pd.DataFrame()
    trans_flux_data = pd.DataFrame()
    for thr in [70]:
        only_trans, only_flux, trans_flux = get_result(thr)
        only_trans_data = pd.concat([only_trans_data, only_trans], axis=0)
        only_flux_data = pd.concat([only_flux_data, only_flux], axis=0)
        trans_flux_data = pd.concat([trans_flux_data, trans_flux], axis=0)
    return only_trans_data, only_flux_data, trans_flux_data


# do things
if __name__ == '__main__':
    print('1')
    n = 70
    msb = pd.read_csv('../data/msbdata.csv')
    flux = pd.read_csv('../output/Glucose_fluxDataset.csv')
    sortedgene, lable = prepare_gene_class(n)
    trans_data = msb.drop(columns=['geneName', 'commonName']).T
    flux_data = flux.drop(columns=['Var1']).T
    ml_trans_data = prepare_mldata(sortedgene, trans_data)
    ml_flux_data = prepare_mldata(sortedgene, flux_data)
    ft = pd.concat([ml_trans_data, ml_flux_data], axis=1)
    logging.info('data loaded')
    #run_ml(ml_flux_data, lable, n, 'flux')
    # run_ml(ft, lable, n, 'ft')
    # run_ml(ml_trans_data, lable, n, 'trans')
# report = classification_report(y, predictions)
