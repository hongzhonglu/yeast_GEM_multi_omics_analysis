import pandas as pd
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, recall_score
from sklearn.neural_network import MLPClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.decomposition import KernelPCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
import logging
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import cross_val_score

logging.basicConfig(level=logging.DEBUG,
format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')


def svm_class(ml_data, lable):
    train_data, test_data, train_lable, test_lable = train_test_split(
        ml_data,
        lable,
        random_state=1,
        train_size=0.8,
        test_size=0.2,
    )
    # train_lable, test_lable = train_test_split(
    #     lable,
    #     random_state=1,
    #     train_size=0.8,
    #     test_size=0.2
    # )

    clf = svm.SVC(
        C=2,
        kernel='rbf',
        gamma=10,
        decision_function_shape='ovr',
        probability=True
    )
    clf.fit(train_data, train_lable)
    scores = cross_val_score(clf, train_data, train_lable, cv=5)
    pre_train_pro = clf.predict_proba(train_data)
    pre_test_pro = clf.predict_proba(test_data)
    pre_pro = np.vstack((pre_train_pro, pre_test_pro))

    pre_train = clf.predict(train_data)
    pre_test = clf.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)
    all_acc = 0.8 * train_acc + 0.2 * test_acc
    pre = np.hstack((pre_train, pre_test))
    true_data = np.hstack((train_lable, test_lable))
    recall_t = recall_score(test_lable, pre_test, average="macro")
    fscore_t = f1_score(test_lable, pre_test, average="macro")
    recall = recall_score(true_data, pre, average="macro")
    fscore = f1_score(true_data, pre, average="macro")
    return train_acc, test_acc, all_acc, recall, fscore, pre, true_data, recall_t, fscore_t, pre_pro


def MLP_class(ml_data, lable):
    train_data, test_data = train_test_split(
        ml_data,
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

    clf = MLPClassifier(
        solver='sgd',
        alpha=1e-5,
        learning_rate='invscaling',
        hidden_layer_sizes=(10, 10),
        max_iter=10,
        random_state=1,
    )
    clf.fit(train_data, train_lable)

    pre_train_pro = clf.predict_proba(train_data)
    pre_test_pro = clf.predict_proba(test_data)
    pre_pro = np.vstack((pre_train_pro, pre_test_pro))

    pre_train = clf.predict(train_data)
    pre_test = clf.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)
    all_acc = 0.8 * train_acc + 0.2 * test_acc
    pre = np.hstack((pre_train, pre_test))
    true_data = np.hstack((train_lable, test_lable))
    recall_t = recall_score(test_lable, pre_test, average="macro")
    fscore_t = f1_score(test_lable, pre_test, average="macro")
    recall = recall_score(true_data, pre, average="macro")
    fscore = f1_score(true_data, pre, average="macro")
    return train_acc, test_acc, all_acc, recall, fscore, pre, true_data, recall_t, fscore_t, pre_pro


def bayes_class(ml_data, lable):
    train_data, test_data, train_lable, test_lable = train_test_split(
        ml_data,
        lable,
        random_state=1,
        train_size=0.8,
        test_size=0.2,
    )
    gnb = GaussianNB()
    gnb.fit(train_data, train_lable)

    pre_train_pro = gnb.predict_proba(train_data)
    pre_test_pro = gnb.predict_proba(test_data)
    pre_pro = np.vstack((pre_train_pro, pre_test_pro))

    pre_train = gnb.predict(train_data)
    pre_test = gnb.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)
    all_acc = 0.8 * train_acc + 0.2 * test_acc
    pre = np.hstack((pre_train, pre_test))
    true_data = np.hstack((train_lable, test_lable))
    recall_t = recall_score(test_lable, pre_test, average="macro")
    fscore_t = f1_score(test_lable, pre_test, average="macro")
    recall = recall_score(true_data, pre, average="macro")
    fscore = f1_score(true_data, pre, average="macro")
    return train_acc, test_acc, all_acc, recall, fscore, pre, true_data, recall_t, fscore_t, pre_pro


def kpca_reduce(ml_data, components):
    transformer = KernelPCA(n_components=components, kernel='rbf')
    X_transformed = transformer.fit_transform(ml_data)
    logging.info('data reduced')
    return X_transformed


def randomf_class(ml_data, lable):
    train_data, test_data, train_lable, test_lable = train_test_split(
        ml_data,
        lable,
        random_state=1,
        train_size=0.8,
        test_size=0.2,
    )
    ranf = RandomForestClassifier()
    ranf.fit(train_data, train_lable)

    pre_train_pro = ranf.predict_proba(train_data)
    pre_test_pro = ranf.predict_proba(test_data)
    pre_pro = np.vstack((pre_train_pro, pre_test_pro))

    pre_train = ranf.predict(train_data)
    pre_test = ranf.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)
    all_acc = 0.8 * train_acc + 0.2 * test_acc
    pre = np.hstack((pre_train, pre_test))
    true_data = np.hstack((train_lable, test_lable))
    recall_t = recall_score(test_lable, pre_test, average="macro")
    fscore_t = f1_score(test_lable, pre_test, average="macro")
    recall = recall_score(true_data, pre, average="macro")
    fscore = f1_score(true_data, pre, average="macro")
    return train_acc, test_acc, all_acc, recall, fscore, pre, true_data, recall_t, fscore_t, pre_pro


def knn_class(ml_data, lable):
    train_data, test_data, train_lable, test_lable = train_test_split(
        ml_data,
        lable,
        random_state=1,
        train_size=0.8,
        test_size=0.2,
    )
    knn = KNeighborsClassifier()
    knn.fit(train_data, train_lable)

    pre_train_pro = knn.predict_proba(train_data)
    pre_test_pro = knn.predict_proba(test_data)
    pre_pro = np.vstack((pre_train_pro, pre_test_pro))

    pre_train = knn.predict(train_data)
    pre_test = knn.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)
    all_acc = 0.8 * train_acc + 0.2 * test_acc
    pre = np.hstack((pre_train, pre_test))
    true_data = np.hstack((train_lable, test_lable))
    recall_t = recall_score(test_lable, pre_test, average="macro")
    fscore_t = f1_score(test_lable, pre_test, average="macro")
    recall = recall_score(true_data, pre, average="macro")
    fscore = f1_score(true_data, pre, average="macro")
    return train_acc, test_acc, all_acc, recall, fscore, pre, true_data, recall_t, fscore_t, pre_pro


def tree_class(ml_data, lable):
    train_data, test_data, train_lable, test_lable = train_test_split(
        ml_data,
        lable,
        random_state=1,
        train_size=0.8,
        test_size=0.2,
    )
    clf = DecisionTreeClassifier()
    clf.fit(train_data, train_lable)

    pre_train_pro = clf.predict_proba(train_data)
    pre_test_pro = clf.predict_proba(test_data)
    pre_pro = np.vstack((pre_train_pro, pre_test_pro))

    pre_train = clf.predict(train_data)
    pre_test = clf.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)
    all_acc = 0.8 * train_acc + 0.2 * test_acc
    pre = np.hstack((pre_train, pre_test))
    true_data = np.hstack((train_lable, test_lable))
    recall_t = recall_score(test_lable, pre_test, average="macro")
    fscore_t = f1_score(test_lable, pre_test, average="macro")
    recall = recall_score(true_data, pre, average="macro")
    fscore = f1_score(true_data, pre, average="macro")
    return train_acc, test_acc, all_acc, recall, fscore, pre, true_data, recall_t, fscore_t, pre_pro


def gradient_class(ml_data, lable):
    train_data, test_data, train_lable, test_lable = train_test_split(
        ml_data,
        lable,
        random_state=1,
        train_size=0.8,
        test_size=0.2,
    )
    clf = GradientBoostingClassifier(learning_rate=0.01, max_depth=8, n_estimators=50)
    clf.fit(train_data, train_lable)

    pre_train_pro = clf.predict_proba(train_data)
    pre_test_pro = clf.predict_proba(test_data)
    pre_pro = np.vstack((pre_train_pro, pre_test_pro))

    pre_train = clf.predict(train_data)
    pre_test = clf.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)
    all_acc = 0.8 * train_acc + 0.2 * test_acc
    pre = np.hstack((pre_train, pre_test))
    true_data = np.hstack((train_lable, test_lable))
    recall_t = recall_score(test_lable, pre_test, average="macro")
    fscore_t = f1_score(test_lable, pre_test, average="macro")
    recall = recall_score(true_data, pre, average="macro")
    fscore = f1_score(true_data, pre, average="macro")
    return train_acc, test_acc, all_acc, recall, fscore, pre, true_data, recall_t, fscore_t, pre_pro


def logistic_class(ml_data, lable):
    train_data, test_data, train_lable, test_lable = train_test_split(
        ml_data,
        lable,
        random_state=1,
        train_size=0.8,
        test_size=0.2,
    )
    clf = LogisticRegression(
        multi_class="multinomial",
        solver="newton-cg",
        max_iter=1000,
    )
    clf.fit(train_data, train_lable)

    pre_train_pro = clf.predict_proba(train_data)
    pre_test_pro = clf.predict_proba(test_data)
    pre_pro = np.vstack((pre_train_pro, pre_test_pro))

    pre_train = clf.predict(train_data)
    pre_test = clf.predict(test_data)
    train_acc = accuracy_score(train_lable, pre_train)
    test_acc = accuracy_score(test_lable, pre_test)
    all_acc = 0.8 * train_acc + 0.2 * test_acc
    pre = np.hstack((pre_train, pre_test))
    true_data = np.hstack((train_lable, test_lable))
    recall_t = recall_score(test_lable, pre_test, average="macro")
    fscore_t = f1_score(test_lable, pre_test, average="macro")
    recall = recall_score(true_data, pre, average="macro")
    fscore = f1_score(true_data, pre, average="macro")
    return train_acc, test_acc, all_acc, recall, fscore, pre, true_data, recall_t, fscore_t, pre_pro


def prepare_mldata(sortedgene, flux_or_trans):
    ml_data = pd.DataFrame(index=sortedgene.index, columns=flux_or_trans.columns)
    for g in sortedgene.index:
        ml_data.loc[g, :] = flux_or_trans.loc[g, :]
    ml_data.index = range(len(ml_data.index))
    return ml_data


def roc(ture_data, pre_pro, n, method, reduction, datatype):
    if n == 50:
        n_class = 9
    elif n == 30:
        n_class = 13
    elif n == 90:
        n_class = 4
    elif n == 70:
        n_class = 5
    lable = label_binarize(ture_data, np.arange(n_class))  # 装换成类似二进制的编码
    # y_score = label_binarize(pre, np.arange(n_class))
    # 1、调用函数计算micro类型的AUC
    # 2、手动计算micro类型的AUC
    # 首先将矩阵y_one_hot和y_score展开，然后计算假正例率FPR和真正例率TPR
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
    plt.savefig('../output/new_roc/ROC and AUC of {} {} {} {}.png'.format(method, n, reduction, datatype),
                dpi=600,
                transparent=True)
    plt.show()
    logging.info('figure saved')


def run_ml(ml_data, lable, n, components, reduction=False):
    svmtrain_acc, svmtest_acc, svmall_acc, svmrecall, svmfscore, svmpre, svmtrue_data, svmrecall_t, svmfscore_t, svmprepro = svm_class(ml_data, lable)
    logging.info('svm finished')
    MLPtrain_acc, MLPtest_acc, MLPall_acc, MLPrecall, MLPfscore, MLPpre, MLPtrue_data, MLPrecall_t, MLPfscore_t, MLPprepro = MLP_class(ml_data, lable)
    logging.info('MLP finished')
    bayestrain_acc, bayestest_acc, bayesall_acc, bayesrecall, bayesfscore, bayespre, bayestrue_data, bayesrecall_t, bayesfscore_t, bayesprepro = bayes_class(ml_data, lable)
    logging.info('bayes finished')
    ranftrain_acc, ranftest_acc, ranfall_acc, ranfrecall, ranffscore, ranfpre, ranftrue_data, ranfrecall_t, ranffscore_t, ranfprepro = randomf_class(ml_data, lable)
    logging.info('randf finished')
    knntrain_acc, knntest_acc, knnall_acc, knnrecall, knnfscore, knnpre, knntrue_data, knnrecall_t, knnfscore_t, knnfprepro = knn_class(ml_data, lable)
    logging.info('knn finished')
    treetrain_acc, treetest_acc, treeall_acc, treerecall, treefscore, treepre, treetrue_data, treerecall_t, treefscore_t, treeprepro = tree_class(ml_data, lable)
    logging.info('tree finished')
    # gratrain_acc, gratest_acc, graall_acc, grarecall = gradient_class(ml_data, lable)
    # logging.info('gra finished')
    logistictrain_acc, logistictest_acc, logisticall_acc, logisticrecall, logisfscore, logispre, logistrue_data, logisrecall_t, logisfscore_t, logisprepro = logistic_class(ml_data, lable)
    logging.info('log finished')
    logging.info('finish computing')
    output = pd.DataFrame(index=['svm_{}'.format(n), 'MLP_{}'.format(n), 'bayes_{}'.format(n),
                                 'randomforest_{}'.format(n), 'knn_{}'.format(n), 'decisiontree_{}'.format(n), 'logistic_{}'.format(n)],
                          columns=['train acc', 'test acc', 'all acc', 'recall', 'fscore', 'pre', 'true_data', 'recall_t', 'fscore_t', 'reduced', 'pre_pro'])
    # output = pd.DataFrame(index=['svm_{}'.format(n), 'bayes_{}'.format(n),
    #                              'randomforest_{}'.format(n), 'knn_{}'.format(n), 'decisiontree_{}'.format(n), 'logistic_{}'.format(n)],
    #                       columns=['train acc', 'test acc', 'all acc', 'recall', 'fscore', 'pre', 'true_data', 'recall_t', 'fscore_t', 'reduced', 'pre_pro'])

    if reduction:
        output.loc['svm_{}'.format(n), :] = [svmtrain_acc, svmtest_acc, svmall_acc, svmrecall, svmfscore, svmpre, svmtrue_data, svmrecall_t, svmfscore_t, 'kpca_{}'.format(components), svmprepro]
        output.loc['MLP_{}'.format(n), :] = [MLPtrain_acc, MLPtest_acc, MLPall_acc, MLPrecall, MLPfscore, MLPpre, MLPtrue_data, MLPrecall_t, MLPfscore_t, 'kpca_{}'.format(components), MLPprepro]
        output.loc['bayes_{}'.format(n), :] = [bayestrain_acc, bayestest_acc, bayesall_acc, bayesrecall, bayesfscore, bayespre, bayestrue_data, bayesrecall_t, bayesfscore_t, 'kpca_{}'.format(components), bayesprepro]
        output.loc['randomforest_{}'.format(n), :] = [ranftrain_acc, ranftest_acc, ranfall_acc, ranfrecall, ranffscore, ranfpre, ranftrue_data, ranfrecall_t, ranffscore_t, 'kpca_{}'.format(components), ranfprepro]
        output.loc['knn_{}'.format(n), :] = [knntrain_acc, knntest_acc, knnall_acc, knnrecall, knnfscore, knnpre, knntrue_data, knnrecall_t, knnfscore_t, 'kpca_{}'.format(components), knnfprepro]
        output.loc['decisiontree_{}'.format(n), :] = [treetrain_acc, treetest_acc, treeall_acc, treerecall, treefscore, treepre, treetrue_data, treerecall_t, treefscore_t, 'kpca_{}'.format(components), treeprepro]
        output.loc['logistic_{}'.format(n), :] = [logistictrain_acc, logistictest_acc, logisticall_acc, logisticrecall, logisfscore, logispre, logistrue_data, logisrecall_t, logisfscore_t, 'kpca_{}'.format(components), logisprepro]

    else:
        output.loc['svm_{}'.format(n), :] = [svmtrain_acc, svmtest_acc, svmall_acc, svmrecall, svmfscore, svmpre, svmtrue_data, svmrecall_t, svmfscore_t, 'off', svmprepro]
        output.loc['MLP_{}'.format(n), :] = [MLPtrain_acc, MLPtest_acc, MLPall_acc, MLPrecall, MLPfscore, MLPpre, MLPtrue_data, MLPrecall_t, MLPfscore_t, 'off', MLPprepro]
        output.loc['bayes_{}'.format(n), :] = [bayestrain_acc, bayestest_acc, bayesall_acc, bayesrecall, bayesfscore, bayespre, bayestrue_data, bayesrecall_t, bayesfscore_t, 'off', bayesprepro]
        output.loc['randomforest_{}'.format(n), :] = [ranftrain_acc, ranftest_acc, ranfall_acc, ranfrecall, ranffscore, ranfpre, ranftrue_data, ranfrecall_t, ranffscore_t, 'off', ranfprepro]
        output.loc['knn_{}'.format(n), :] = [knntrain_acc, knntest_acc, knnall_acc, knnrecall, knnfscore, knnpre, knntrue_data, knnrecall_t, knnfscore_t, 'off', knnfprepro]
        output.loc['decisiontree_{}'.format(n), :] = [treetrain_acc, treetest_acc, treeall_acc, treerecall, treefscore, treepre, treetrue_data, treerecall_t, treefscore_t, 'off', treeprepro]
        output.loc['logistic_{}'.format(n), :] = [logistictrain_acc, logistictest_acc, logisticall_acc, logisticrecall, logisfscore, logispre, logistrue_data, logisrecall_t, logisfscore_t, 'off', logisprepro]
    return output


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
                + [2] * 148 + [3] * 92\
                + [4] * 110
    return sortedgene, lable


def get_result(n, components=0):
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

    if components:
        # only with trans
        ml_trans_data_re = kpca_reduce(ml_trans_data, components)
        only_trans = run_ml(ml_trans_data_re, lable, n, components, reduction=True)
        for m in only_trans.index.tolist():
            roc(only_trans.loc[m, 'true_data'], only_trans.loc[m, 'pre_pro'], n, m, only_trans.loc[m, 'reduced'], 'only_trans')
        logging.info('trans finished')
        # only with flux
        ml_flux_data_re = kpca_reduce(ml_flux_data, components)
        only_flux = run_ml(ml_flux_data_re, lable, n, components, reduction=True)
        for m in only_flux.index.tolist():
            roc(only_flux.loc[m, 'true_data'], only_flux.loc[m, 'pre_pro'], n, m, only_flux.loc[m, 'reduced'], 'only_flux')
        logging.info('flux finished')
        # flux and trans
        ft = kpca_reduce(ft, components)
        trans_flux = run_ml(ft, lable, n, components, reduction=True)
        for m in only_flux.index.tolist():
            roc(trans_flux.loc[m, 'true_data'], trans_flux.loc[m, 'pre_pro'], n, m, trans_flux.loc[m, 'reduced'], 'trans_flux')
        logging.info('all finished')

    else:
        # only with trans
        only_trans = run_ml(ml_trans_data, lable, n, components)
        for m in only_trans.index.tolist():
            roc(only_trans.loc[m, 'true_data'], only_trans.loc[m, 'pre_pro'], n, m, only_trans.loc[m, 'reduced'], 'only_trans')
        logging.info('trans finished')
        # only with flux
        only_flux = run_ml(ml_flux_data, lable, n, components)
        for m in only_flux.index.tolist():
            roc(only_flux.loc[m, 'true_data'], only_flux.loc[m, 'pre_pro'], n, m, only_flux.loc[m, 'reduced'], 'only_flux')
        logging.info('flux finished')
        # flux and trans
        trans_flux = run_ml(ft, lable, n, components)
        for m in only_flux.index.tolist():
            roc(trans_flux.loc[m, 'true_data'], trans_flux.loc[m, 'pre_pro'], n, m, trans_flux.loc[m, 'reduced'], 'trans_flux')
        logging.info('all finished')
    return only_trans, only_flux, trans_flux


def from_result_to_output(components):
    only_trans_data = pd.DataFrame()
    only_flux_data = pd.DataFrame()
    trans_flux_data = pd.DataFrame()
    for thr in [30, 50, 70, 90]:
        only_trans, only_flux, trans_flux = get_result(thr, components)
        only_trans_data = pd.concat([only_trans_data, only_trans], axis=0)
        only_flux_data = pd.concat([only_flux_data, only_flux], axis=0)
        trans_flux_data = pd.concat([trans_flux_data, trans_flux], axis=0)
    return only_trans_data, only_flux_data, trans_flux_data


# do things
if __name__ == '__main__':
    for rd in list(range(5)):
        print(rd)
        only_trans_data_off, only_flux_data_off, trans_flux_data_off = from_result_to_output(components=0)
        only_trans_data_200, only_flux_data_200, trans_flux_data_200 = from_result_to_output(components=200)
        only_trans_data_500, only_flux_data_500, trans_flux_data_500 = from_result_to_output(components=500)
        only_trans_data_800, only_flux_data_800, trans_flux_data_800 = from_result_to_output(components=800)
        only_trans_data_1000, only_flux_data_1000, trans_flux_data_1000 = from_result_to_output(components=1000)

        trans = pd.concat([only_trans_data_off, only_trans_data_200,
                           only_trans_data_500, only_trans_data_800,
                           only_trans_data_1000], axis=0)

        flux = pd.concat([only_flux_data_off, only_flux_data_200,
                          only_flux_data_500, only_flux_data_800,
                          only_flux_data_1000], axis=0)

        trans_flux = pd.concat([trans_flux_data_off, trans_flux_data_200,
                                trans_flux_data_500, trans_flux_data_800,
                                trans_flux_data_1000], axis=0)

        writer = pd.ExcelWriter('../output/ML_result{}.xlsx'.format(rd))
        sheet_name1 = 'only trans'
        trans.to_excel(
            excel_writer=writer,
            sheet_name=sheet_name1)
        sheet_name2 = 'only flux'
        flux.to_excel(
            excel_writer=writer,
            sheet_name=sheet_name2)
        sheet_name3 = 'trans and flux'
        trans_flux.to_excel(
            excel_writer=writer,
            sheet_name=sheet_name3)
        writer.save()
