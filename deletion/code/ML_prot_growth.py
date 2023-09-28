import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import re
import time

def get_col(orf, data):
    for col in data.columns:
        if len(re.compile('{}_'.format(orf)).findall(col)) != 0:
            break
    return col


def train(data, sc):
    X_train, X_test, y_train, y_test = train_test_split(data.T,
                                                        sc,
                                                        test_size=0.2,
                                                        random_state=0)

    from sklearn.preprocessing import StandardScaler

    sc = StandardScaler()
    X_train = sc.fit_transform(X_train)
    X_test = sc.transform(X_test)

    from sklearn.ensemble import RandomForestRegressor

    regressor = RandomForestRegressor(n_estimators=200, random_state=0)
    regressor.fit(X_train, y_train)
    y_pred = regressor.predict(X_test)

    from sklearn import metrics
    from sklearn.metrics import r2_score
    r2_score(y_test, y_pred)

    print('Mean Absolute Error:', metrics.mean_absolute_error(y_test, y_pred))
    print('Mean Squared Error:', metrics.mean_squared_error(y_test, y_pred))
    print('Root Mean Squared Error:',
          np.sqrt(metrics.mean_squared_error(y_test, y_pred)))
    print('r2:',
          r2_score(y_test, y_pred))


def get_train(target, data):
    common = [p for p in target if p in list(prot.index)]
    tar_prot = pd.DataFrame(index=common, columns=data.columns)
    for t in common:
        tar_prot.loc[t, :] = data.loc[t, :]
    return tar_prot


prot = pd.read_csv('../data/yeast5k_noimpute_wide_without_target.csv')
growth = pd.read_csv('../data/yeast5k_growthrates_byORF_without_target.csv', encoding='latin-1')
target1 = ['YNL069C', 'YPL198W', 'YPL217C', 'YIL133C', 'YGL123W', 'YNL096C', 'YJL191W', 'YMR121C', 'YOR145C', 'YPL043W', 'YDL153C', 'YDR365C', 'YJR002W', 'YDR312W', 'YKL009W', 'YLR009W', 'YLR051C']
target2 = ['YML063W','YMR143W','YPL081W','YKL081W','YDL130W','YPL274W','YOR063W','YBL087C','YER091C','YDR194C','YLR441C','YNL069C','YBR031W','YLR325C','YBR189W','YNL178W','YPR041W','YBR181C','YBR191W','YPL198W','YGL031C','YGL076C','YML024W','YPL131W','YGR204W','YGL030W','YPL160W','YPL143W','YGL135W','YHL015W','YPL217C','YOR234C','YDL059C','YHL033C','YIR034C','YHR020W','YDR408C','YIL133C','YBR122C','YGL123W','YLR364W','YFL036W','YOL120C','YGR214W','YHR010W','YNL096C','YLR029C','YGR034W','YJL191W','YLR075W','YNL067W','YDL202W','YDL191W','YMR121C','YMR260C','YAL035W','YOR361C','YHR024C','YLR264W','YBR282W','YJR123W','YIL052C','YBR267W','YLR312W-A','YLR363W-A','YGL200C','YOL121C','YLR344W','YOR254C','YGR118W','YKR057W','YPR033C','YGL236C','YML026C','YLL045C','YFR030W','YLR439W','YLR249W','YHR068W','YKL001C','YJR045C','YCR046C','YMR157C','YPL249C-A','YDR405W','YDR232W','YIR030C','YDR175C','YAL036C','YNL306W','YLR448W','YOR145C','YOR091W','YOL039W','YOR206W','YGR094W','YPL266W','YJR137C','YKL170W','YLR340W','YGL189C','YPR100W','YOR096W','YGR174C','YOL040C','YPL043W','YMR024W','YHL004W','YDL153C','YNL209W','YHR208W','YOR065W','YGL026C','YOR187W','YDR321W','YBR084C-A','YDR365C','YOR204W','YAL012W','YJL104W','YNR050C','YJR002W','YMR203W','YDR462W','YOL097C','YLR333C','YDR312W','YCL030C','YDR382W','YDR064W','YKL009W','YMR142C','YNL124W','YBR248C','YDL166C','YEL054C','YHR203C','YML073C','YLR009W','YGR191W','YOR182C','YER165W','YLR051C','YBR262C','YFR031C-A','YBR143C','YHR019C']
# prot = pd.read_csv('../data/yeast5k_noimpute_wide.csv')
# growth = pd.read_csv('../data/yeast5k_growthrates_byORF.csv', encoding='latin-1')
sc = growth.loc[:, 'SC']
growth.index = growth.loc[:, 'orf']
prot.index = prot.loc[:, 'sys_name']
data = pd.DataFrame(index=prot.index, columns=[get_col(i, prot) for i in growth.index])
c = 0

for i in growth.index:
    print(c)
    col = get_col(i, prot)
    data.loc[:, col] = prot.loc[:, col]
    c += 1


tar_prot1 = get_train(target1, data)
tar_prot2 = get_train(target2, data)
print('target1')
st = time.time()
train(tar_prot1, sc)
ed = time.time()
print('Training time:{}'.format((ed-st)/3600))
print('target2')
st = time.time()
train(tar_prot2, sc)
ed = time.time()
print('Training time:{}'.format((ed-st)/3600))
