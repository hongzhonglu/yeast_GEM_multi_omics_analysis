# find which rxn is closed by GIMME


import cobra
import logging
from tqdm import trange
import pandas as pd


def get_closed_rxn(y9, scmodel):
    closed_rxn = [r.id for r in y9.reactions if r not in scmodel.reactions]
    return closed_rxn


def read_model(filepath):
    try:
        scmodel = cobra.io.read_sbml_model(filepath)
    except OSError:
        scmodel = ''
    return scmodel


def count_rxn(y9, condition):
    dic = {}
    for r in y9.reactions:
        count = 0
        for c in condition.values():
            if r.id in c:
                count += 1
        dic[r.id] = 100 * count/len(condition)
    return dic


# do things
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')

    y9 = cobra.io.read_sbml_model('../data/yeast-GEM.xml')
    stress = {}
    unstress = {}
    for sp in trange(1, 95):
        # stress
        filepath_s = '../output/scmodel/BY4741_Stressed_1784042160_C50_BY{}.xml'.format(str(sp))
        scmodel_s = read_model(filepath_s)
        if scmodel_s != '':
            closed_rxn_s = get_closed_rxn(y9, scmodel_s)
            stress['Stressed_BY{}'.format(str(sp))] = closed_rxn_s
        # unstress
        filepath_u = '../output/scmodel/BY4741_Unstressed_1784042189_BY{}.xml'.format(str(sp))
        scmodel_u = read_model(filepath_u)
        if scmodel_u != '':
            closed_rxn_u = get_closed_rxn(y9, scmodel_u)
            unstress['Unstressed_BY{}'.format(str(sp))] = closed_rxn_u

    # count
    count_s = count_rxn(y9, stress)
    count_u = count_rxn(y9, unstress)
    sub = pd.read_table('../data/uniqueSubsystems.tsv')
    y9excel = pd.read_excel(r"C:\Users\yuhuzhouye\Desktop\yeast9_w\yeast-GEM.xlsx")
    y9excel.index = y9excel.loc[:, 'ID']
    sub.index = sub.loc[:, 'ID']
    diff = pd.DataFrame(columns=['stress', 'unstress', 'subsystem', 'Equation'])
    for rr in count_u.keys():
        if (max(count_s[rr], count_u[rr]) >= 2 * min(count_s[rr], count_u[rr])) and \
           (max(count_s[rr], count_u[rr]) > 50):
            diff.loc[rr, 'stress'] = count_s[rr]
            diff.loc[rr, 'unstress'] = count_u[rr]
            diff.loc[rr, 'subsystem'] = sub.loc[rr, 'subsystem_unique_@2021_12']
            diff.loc[rr, 'Equation'] = y9excel.loc[rr, 'Equation']
    diff.to_excel('../output/diff.xlsx')