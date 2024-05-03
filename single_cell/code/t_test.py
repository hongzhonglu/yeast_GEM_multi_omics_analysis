##########################################
# get "t-test.jpg"
# compute the p value of the same reaction flux between stress and unstress condition
# first run "get_pfba_flux.py" to compute flux data "pfba_flux.xlsx"

from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


pathway = pd.read_table('../data/uniqueSubsystems.tsv')
pathway.index = pathway.loc[:, 'ID']
t_data = pd.read_excel('../output/pfba_flux.xlsx')

t_data = t_data.iloc[:, 3:]

co = t_data.columns.tolist()
for z in co:
    if sum(t_data.loc[:, z]) == 0:
        t_data = t_data.drop(z, axis=1)

stress = t_data.loc[0:79, :]
unstress = t_data.loc[80:, :]
none_zero_l = stress.columns.tolist()
rxn_p = {}
p_xiao = {}
p_da = {}
for te in none_zero_l:
    # 独立2个样本t检验
    sample1 = np.asarray(stress.loc[:, te])
    sample2 = np.asarray(unstress.loc[:, te])
    r = stats.ttest_ind(sample1, sample2)
    rxn_p[te] = r.pvalue
    if r.pvalue <= 0.05:
        p_xiao[te] = r.pvalue
    if r.pvalue > 0.05:
        p_da[te] = r.pvalue

stress_mean = stress.mean()
unstress_mean = unstress.mean()
#
# mean flux of each rxn(salt stress)
stress_mean.to_excel('../output/stress_flux_mean.xlsx')
# mean flux of each rxn(unstress)
unstress_mean.to_excel('../output/unstress_flux_mean.xlsx')
# # rxn(p < 0.05) subsystem
# xiao_path_S.to_excel('../output/path_analysis.xlsx')
# fenlei.to_excel('../output/t-test_0.05.xlsx')

# figure
name_l = ['p value > 0.05', 'p value ≤ 0.05']
values = [len(none_zero_l) - len(p_xiao), len(p_xiao)]
font1 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 23,
         }
font2 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 20,
         }
font3 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 12,
         }

plt.subplots(figsize=(8, 6))
plt.pie(values,
        labeldistance=1,
        wedgeprops={'linewidth': 3, 'edgecolor': 'white'},
        colors=['#64F0E1', '#F09165'],
        textprops=font2)
plt.title('t-test of nonzero flux', fontdict=font1)
plt.legend(name_l,
           bbox_to_anchor=(0.9, 0.9),
           prop=font3)
plt.savefig('../output/t-test.tif')
plt.show()
