from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

pathway = pd.read_table(r'singlecell/data/Rxn_unique_subsystem.tsv')
pathway.index = pathway.loc[:, 'ID']
t_data = pd.read_excel('singlecell/output/pfba_flux.xlsx')
t_data = t_data.drop(0, axis=1)
t_data = t_data.drop('metnum/20', axis=1)
t_data = t_data.drop('rxnnum/20', axis=1)

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
    r = stats.ttest_ind(sample1, sample2, )
    rxn_p[te] = r.pvalue
    if r.pvalue <= 0.05:
        p_xiao[te] = r.pvalue
    if r.pvalue > 0.05:
        p_da[te] = r.pvalue

xiao_path = {k: pathway.loc[k, 'subsystem_unique_@2021_12'] for k in p_xiao.keys()}
xiao_path_l = [i for i in xiao_path.values()]
set01 = set(xiao_path_l)
dict01 = {}
for item in set01:
    dict01.update({item:xiao_path_l.count(item)})
fenlei = pd.Series(dict01)
print('总计：' + str(len(none_zero_l)) + "个非零通量\n"
      + str(len(p_xiao)) + '个差异显著的通量')

xiao_path_S = pd.Series(xiao_path)
stress_mean = stress.mean()
unstress_mean = unstress.mean()

stress_mean.to_excel(r'singlecell/output/stress_flux_mean.xlsx')
unstress_mean.to_excel(r'singlecell/output/unstress_flux_mean.xlsx')
xiao_path_S.to_excel(r'singlecell/output/path_analysis.xlsx')
fenlei.to_excel(r'singlecell/output/t-test_0.05.xlsx')

fenlei_T = pd.read_excel(r'singlecell/output/t-test_0.05.xlsx',
                         sheet_name='Sheet2')
bar_data = fenlei_T.loc[0:12, 0]
bar_name = fenlei_T.loc[0:12, 'subsystem']

plt.rcParams['figure.figsize'] = (8.0, 4.0)
plt.barh(range(len(bar_data)), bar_data, tick_label=bar_name)
plt.show()