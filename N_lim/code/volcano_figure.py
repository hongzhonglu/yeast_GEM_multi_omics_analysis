import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
from scipy import stats

os.chdir('..')
os.chdir('..')
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
p = {}
for te in none_zero_l:
    # 独立2个样本t检验
    sample1 = np.asarray(stress.loc[:, te])
    sample2 = np.asarray(unstress.loc[:, te])
    r = stats.ttest_ind(sample1, sample2)
    if r.pvalue <= 0.05:
        p[te] = r.pvalue
    if r.pvalue > 0.05:
        p[te] = r.pvalue

p_value = [v for v in p.values()]
p_key = [k for k in p.keys()]
stress_mean = stress.mean()
unstress_mean = unstress.mean()

s = np.asarray(stress_mean)
us = np.asarray(unstress_mean)
divi = []
for ss in range(len(s)):
    if s[ss] == 0:
        divi.append(1000000)
    elif us[ss] == 0:
        divi.append(2000000)
    else:
        divi.append(abs(s[ss]/us[ss]))

tempfold = {}
tempp = {}
s0 = {}
us0 = {}
for l in range(len(divi)):
    if divi[l] < 1000000:
        tempfold[p_key[l]] = divi[l]
        tempp[p_key[l]] = p_value[l]
    if divi[l] == 1000000:
        s0[p_key[l]] = 's0'
    if divi[l] == 2000000:
        us0[p_key[l]] = 'us0'

log_p = -np.log10([pp for pp in tempp.values()])
log_foldchange = np.log2([ff for ff in tempfold.values()])

x_threshold=0.58
y_threshold=1.3
result = pd.DataFrame(columns=['x', 'y'])
result.loc[:, 'x'] = log_foldchange
result.loc[:, 'y'] = log_p
result['group'] = 'black'
result.loc[(result.x > x_threshold)&(result.y > y_threshold),'group'] = 'tab:red' #x=-+x_threshold直接截断
result.loc[(result.x < -x_threshold)&(result.y > y_threshold),'group'] = 'tab:blue' #x=-+x_threshold直接截断
result.loc[result.y < y_threshold,'group'] = 'dimgrey' #阈值以下点为灰色
print(result.head())

#确定坐标轴显示范围
xmin=-8
xmax=10
ymin=-1
ymax=10

#绘制散点图
fig = plt.figure(figsize=plt.figaspect(7/6)) #确定fig比例（h/w）
ax = fig.add_subplot()
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title='')
ax.scatter(result['x'], result['y'], s=2, c=result['group'])
# ax.set_ylabel('-Log10(Q value)',fontweight='bold')
# ax.set_xlabel('Log2 (fold change)',fontweight='bold')
ax.spines['right'].set_visible(False) #去掉右边框
ax.spines['top'].set_visible(False) #去掉上边框

#水平和竖直线
ax.vlines(-x_threshold, ymin, ymax, color='dimgrey',linestyle='dashed', linewidth=1) #画竖直线
ax.vlines(x_threshold, ymin, ymax, color='dimgrey',linestyle='dashed', linewidth=1) #画竖直线
ax.hlines(y_threshold, xmin, xmax, color='dimgrey',linestyle='dashed', linewidth=1) #画竖水平线

ax.set_xticks(range(-8,12,2)) #设置x轴刻度起点和步长
ax.set_yticks(range(-1,10,2)) #设置y轴刻度起点和步长
plt.ylabel('-log10(p value)',
           size=12,
           fontproperties='Arial',
           )
plt.xlabel('log2(flux fold changge)',
           size=12,
           fontproperties='Arial'
           )
plt.show()

output = pd.DataFrame(columns=['foldchange', 'p', 'logfold', 'logp'],
                      index=list(tempfold.keys()))

output.loc[:, 'foldchange'] = list(tempfold.values())
output.loc[:, 'p'] = list(tempp.values())
output.loc[:, 'logfold'] = [i for i in log_foldchange]
output.loc[:, 'logp'] = [j for j in log_p]
