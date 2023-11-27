def draw(true_data1, pre_pro1, true_data2, pre_pro2, y, x=1):
    lable1 = label_binarize(true_data1, np.arange(5))  # 装换成类似二进制的编码
    fpr1, tpr1, thresholds1 = metrics.roc_curve(lable1.ravel(), pre_pro1.ravel())
    auc1 = metrics.auc(fpr1, tpr1)

    lable = label_binarize(true_data2, np.arange(5))  # 装换成类似二进制的编码
    fpr, tpr, thresholds = metrics.roc_curve(lable.ravel(), pre_pro2.ravel())
    auc = metrics.auc(fpr, tpr)
    # 绘图
    font3 = {'family': 'Arial',
             'weight': 'normal',
             'size': 7,
             }
    if x:
        plt.plot(fpr1, tpr1, c='#467897', lw=1, alpha=1, label=u'AUC=%.3f' % auc1)
        plt.plot(fpr, tpr, c='#E7CD79', lw=1, alpha=1, label=u'AUC=%.3f' % auc)
    else:
        plt.plot(fpr1, tpr1, c='#467897', lw=1, alpha=1, label=u'Transcriptome AUC=%.3f' % auc1)
        plt.plot(fpr, tpr, c='#E7CD79', lw=1, alpha=1, label=u'Fluxomic AUC=%.3f' % auc)

    plt.plot((0, 1), (0, 1), c='#808080', lw=1, ls='--', alpha=1)
    plt.xlim((-0.01, 1.02))
    plt.ylim((-0.01, 1.02))
    plt.xticks(np.arange(0, 1.1, 0.5), fontproperties='Arial', size=6)
    if y:
        plt.yticks(np.arange(0, 1.1, 0.2), fontproperties='Arial', size=6)
    else:
        plt.yticks([])
    # plt.grid(b=True, ls=':')
    plt.legend(loc='lower right', fancybox=True, framealpha=False, prop=font3)


font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 18,
         }
plt.figure(figsize=(6, 3))
plt.subplot(2, 3, 1)
draw(only_trans_data_500.loc['randomforest_70', 'true_data'], only_trans_data_500.loc['randomforest_70', 'pre_pro'],
     only_flux_data_500.loc['randomforest_70', 'true_data'], only_flux_data_500.loc['randomforest_70', 'pre_pro'], y=1,
     x=0)
plt.subplot(2, 3, 2)
draw(only_trans_data_500.loc['MLP_70', 'true_data'], only_trans_data_500.loc['MLP_70', 'pre_pro'],
     only_flux_data_500.loc['MLP_70', 'true_data'], only_flux_data_500.loc['MLP_70', 'pre_pro'], y=0)
plt.subplot(2, 3, 3)
draw(only_trans_data_500.loc['svm_70', 'true_data'], only_trans_data_500.loc['svm_70', 'pre_pro'],
     only_flux_data_500.loc['svm_70', 'true_data'], only_flux_data_500.loc['svm_70', 'pre_pro'], y=0)
plt.subplot(2, 3, 4)
draw(only_trans_data_500.loc['bayes_70', 'true_data'], only_trans_data_500.loc['bayes_70', 'pre_pro'],
     only_flux_data_500.loc['bayes_70', 'true_data'], only_flux_data_500.loc['bayes_70', 'pre_pro'], y=1)
plt.subplot(2, 3, 5)
draw(only_trans_data_500.loc['knn_70', 'true_data'], only_trans_data_500.loc['knn_70', 'pre_pro'],
     only_flux_data_500.loc['knn_70', 'true_data'], only_flux_data_500.loc['knn_70', 'pre_pro'], y=0)
plt.subplot(2, 3, 6)
draw(only_trans_data_500.loc['logistic_70', 'true_data'], only_trans_data_500.loc['logistic_70', 'pre_pro'],
     only_flux_data_500.loc['logistic_70', 'true_data'], only_flux_data_500.loc['logistic_70', 'pre_pro'], y=0)

plt.tight_layout()
plt.savefig('../output/ROC and AUC.png',
            dpi=600,
            transparent=True)
plt.show()
