def draw(true_data1, pre_pro1, true_data2, pre_pro2, true_data3, pre_pro3, y, x=1):
    lable1 = label_binarize(true_data1, np.arange(5))  # 装换成类似二进制的编码
    fpr1, tpr1, thresholds1 = metrics.roc_curve(lable1.ravel(), pre_pro1.ravel())
    auc1 = metrics.auc(fpr1, tpr1)

    lable = label_binarize(true_data2, np.arange(5))  # 装换成类似二进制的编码
    fpr, tpr, thresholds = metrics.roc_curve(lable.ravel(), pre_pro2.ravel())
    auc = metrics.auc(fpr, tpr)

    lable3 = label_binarize(true_data3, np.arange(5))  # 装换成类似二进制的编码
    fpr3, tpr3, thresholds3 = metrics.roc_curve(lable3.ravel(), pre_pro3.ravel())
    auc3 = metrics.auc(fpr3, tpr3)

    # 绘图
    font3 = {'family': 'Arial',
             'weight': 'normal',
             'size': 7,
             }
    if x:
        plt.plot(fpr1, tpr1, c='#EA8379', lw=3, alpha=1)
        plt.plot(fpr, tpr, c='#7DAEE0', lw=3, alpha=1)
        plt.plot(fpr3, tpr3, c='#B395BD', lw=3, alpha=0.8)
    else:
        plt.plot(fpr1, tpr1, c='#EA8379', lw=1, alpha=1, label=u'AUC=%.3f' % auc1)
        plt.plot(fpr, tpr, c='#7DAEE0', lw=1, alpha=1, label=u'AUC=%.3f' % auc)
        plt.plot(fpr3, tpr3, c='#B395BD', lw=1, alpha=1, label=u'AUC=%.3f' % auc3)
        plt.legend(prop={'weight': 'bold'})
    plt.plot((0, 1), (0, 1), c='#808080', lw=1, ls='--', alpha=1)
    plt.xlim((-0.01, 1.02))
    plt.ylim((-0.01, 1.02))
    plt.xticks(np.arange(0, 1.1, 0.5), fontproperties='Arial', size=20, fontweight='bold')
    if y:
        plt.yticks(np.arange(0, 1.1, 0.5), fontproperties='Arial', size=20, fontweight='bold')
    else:
        plt.yticks([])
    # plt.grid(b=True, ls=':')
    #plt.legend(loc='lower right', fancybox=True, framealpha=False, prop=font3)


font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 14,
         }
font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 16,
         }
fig, ax = plt.subplots(figsize=(8, 6))# plt.figure(figsize=(6, 3))

# plt.subplot(2, 3, 1)
# draw(true_data_rf_t, pre_pro_rf_t,
#      true_data_rf_f, pre_pro_rf_f,
#      true_data_rf_ft, pre_pro_rf_ft,
#      y=1, x=0)
# plt.subplot(2, 3, 2)
# draw(true_data_mlp_t, pre_pro_mlp_t,
#      true_data_mlp_f, pre_pro_mlp_f,
#      true_data_mlp_ft, pre_pro_mlp_ft,
#      y=0)
# plt.subplot(2, 3, 3)
# draw(true_data_svm_t, pre_pro_svm_t,
#      true_data_svm_f, pre_pro_svm_f,
#      true_data_svm_ft, pre_pro_svm_ft,
#      y=0)
# plt.subplot(2, 3, 4)
draw(true_data_gnb_t, pre_pro_gnb_t,
     true_data_gnb_f, pre_pro_gnb_f,
     true_data_gnb_ft, pre_pro_gnb_ft,
     y=1, x=1)
# plt.subplot(2, 3, 5)
# draw(true_data_knn_t, pre_pro_knn_t,
#      true_data_knn_f, pre_pro_knn_f,
#      true_data_knn_ft, pre_pro_knn_ft,
#      y=0)
# plt.subplot(2, 3, 6)
# draw(true_data_log_t, pre_pro_log_t,
#      true_data_log_f, pre_pro_log_f,
#      true_data_log_ft, pre_pro_log_ft,
#      y=0)
# plt.xlabel('False Positive Rate', fontdict=font1)
# plt.ylabel('True Positive Rate', fontdict=font1)
# plt.grid(b=True, ls=':')
ax.spines['top'].set_linewidth(2)    # 顶部边框
ax.spines['right'].set_linewidth(2)  # 右侧边框
ax.spines['bottom'].set_linewidth(2) # 底部边框
ax.spines['left'].set_linewidth(2)   # 左侧边框

plt.tight_layout()
# plt.savefig('./output/ROC and AUC.png',
#             dpi=600,
#             transparent=True)
plt.savefig('../output/ROC and AUC.jpg',
            dpi=600)

plt.show()
