import pandas as pd
import cobra
from cobra.io import read_sbml_model
import numpy as np
import re
import time
import matplotlib.pyplot as plt
from cobra.flux_analysis import pfba


def normalization(data):
    # return ((1 - 0.0001) * (data - np.min(data))) / (np.max(data) - np.min(data)) + 0.0001
    return (data - np.min(data)) / (np.max(data) - np.min(data))
    # return (data - np.mean(data)) / (np.std(data))


def process_data(data, prot):
    # def sigmoid(x):
    #     s = 1 / (1 + np.exp(-x))
    #     return s
    for col in data.columns:
        data.loc[:, col] = normalization(np.log2(prot.loc[:, col]))
        # data.loc[:, col] = sigmoid(prot.loc[:, col])
    return data


def get_prot_conc(tempgene, sgd, protdict, condition):
    gene1 = re.compile(r'Y.{6}').findall(tempgene)
    gene2 = re.compile(r'Y.{6}-\w').findall(tempgene)
    gene3 = re.compile(r'Q.{4}').findall(tempgene)
    gene = gene1 + gene2 + gene3
    templist = []
    for t in gene:
        try:
            t = sgd.loc[t, 'Standard_name']
            templist.append(protdict[condition].loc[t, 'prot'])
        except KeyError or ValueError:
            pass
    return templist


def oneone(reaction, prot_data):
    try:
        # reaction.lower_bound = reaction.lower_bound * prot_data.loc[reaction.gene_reaction_rule]
        # reaction.upper_bound = reaction.upper_bound * prot_data.loc[reaction.gene_reaction_rule]
        reaction.lower_bound = -np.power(abs(reaction.lower_bound), prot_data.loc[reaction.gene_reaction_rule])
        reaction.upper_bound = np.power(reaction.upper_bound, prot_data.loc[reaction.gene_reaction_rule])

    except KeyError:
        pass


def get_gene(gpr):
    gene1 = re.compile(r'Y.{6}').findall(gpr)
    gene2 = re.compile(r'Y.{6}-\w').findall(gpr)
    gene3 = re.compile(r'Q.{4}').findall(gpr)
    gene = gene1 + gene2 + gene3
    return gene


def get_constrain(gpr, prot_data):
    gene = get_gene(gpr)
    val = []
    for g in gene:
        try:
            val.append(prot_data.loc[g])
        except KeyError:
            val.append(1)
    return val


def est_and_or(gpr):
    re_and = re.compile('and')
    re_or = re.compile('or')
    com_and = re_and.findall(gpr)
    com_or = re_or.findall(gpr)
    if len(com_or) == 0 and len(com_and) == 0:
        relation = 'oneone'
    if len(com_or) == 0 and len(com_and) != 0:
        relation = 'complex'
    if len(com_or) != 0 and len(com_and) == 0:
        relation = 'isoenzymes'
    if len(com_or) != 0 and len(com_and) != 0:
        relation = 'complex and isoenzyme'
    return relation


def constrain_bounds(tmpmodel, prot_data):
    for r in tmpmodel.reactions:
        relation = est_and_or(r.gene_reaction_rule)
        if relation == 'oneone':
            oneone(r, prot_data)
        if relation == 'complex':
            val = get_constrain(r.gene_reaction_rule, prot_data)
            # r.lower_bound = r.lower_bound * np.min(val)
            # r.upper_bound = r.upper_bound * np.min(val)
            r.lower_bound = -np.power(abs(r.lower_bound), np.min(val))
            r.upper_bound = np.power(r.upper_bound, np.min(val))
        if relation == 'isoenzymes':
            val = get_constrain(r.gene_reaction_rule, prot_data)
            # r.lower_bound = r.lower_bound * np.max(val)
            # r.upper_bound = r.upper_bound * np.max(val)
            r.lower_bound = -np.power(abs(r.lower_bound), np.sum(val))
            r.upper_bound = np.power(r.upper_bound, np.sum(val))
        if relation == 'complex and isoenzyme':
            com_iso = r.gene_reaction_rule.split('or')
            tmpval = []
            for com in com_iso:
                tmpval.append(np.min(get_constrain(com, prot_data)))
            val = np.max(tmpval)
            # r.lower_bound = r.lower_bound * np.max(val)
            # r.upper_bound = r.upper_bound * np.max(val)
            r.lower_bound = -np.power(abs(r.lower_bound), np.sum(val))
            r.upper_bound = np.power(r.upper_bound, np.sum(val))
    return tmpmodel


def compute_u(model, prot_data):
    tmpmodel = model.copy()
    tmpmodel = change_sc_medium(tmpmodel)
    tmpmodel = constrain_bounds(tmpmodel, prot_data)
    u = tmpmodel.optimize().objective_value
    # u = pfba(tmpmodel).fluxes.loc['r_2111']
    return u


def change_sc_medium(model):
    # Biotin
    model.reactions.get_by_id('r_1671').lower_bound = -1
    # CoA
    model.reactions.get_by_id('r_1548').lower_bound = -1
    # Inositol
    model.reactions.get_by_id('r_1947').lower_bound = -1
    # nicotinate
    model.reactions.get_by_id('r_1967').lower_bound = -1
    # aminobenzoate
    model.reactions.get_by_id('r_1604').lower_bound = -1
    # pyridoxine
    model.reactions.get_by_id('r_2028').lower_bound = -1
    # riboflavin
    model.reactions.get_by_id('r_2038').lower_bound = -1
    # Thiamine
    model.reactions.get_by_id('r_2067').lower_bound = -1
    medium = model.medium
    medium['r_1671'] = 1.0
    medium['r_1548'] = 1.0
    medium['r_1947'] = 1.0
    medium['r_1947'] = 1.0
    medium['r_1967'] = 1.0
    medium['r_1604'] = 1.0
    medium['r_2028'] = 1.0
    medium['r_2038'] = 1.0
    medium['r_2067'] = 1.0
    medium['r_1714'] = 10.0
    model.medium = medium
    return model


def get_col(orf, data):
    for col in data.columns:
        if len(re.compile(orf, re.IGNORECASE).findall(col)) != 0:
            break
    return col


if __name__ == '__main__':
    model = read_sbml_model(r"C:\Users\yuhuzhouye\Desktop\yeast9_w\yesat9_temp\yeast-GEM.xml")
    prot = pd.read_csv('../data/yeast5k_noimpute_wide.csv')
    exp = pd.read_csv('../data/yeast5k_growthrates_byORF.csv', encoding='latin-1')
    prot.index = prot.iloc[:, 0]

    data = pd.DataFrame(index=prot.index, columns=list(prot.columns)[2: -1])
    data = process_data(data, prot)
    growth = []
    exp_growth = []
    ct = 0
    # for i in range(len(data.columns)):
    start = time.time()
    exp.index = exp.loc[:, 'orf']
    for i in exp[:]['orf']:
        col = get_col(i, data)
        prot_data = data.loc[:, col]
        u = compute_u(model, prot_data)
        growth.append(u)
        ct += 1
        print(ct)
        print(u)
        exp_growth.append(exp.loc[i, 'SC'])
        print(exp.loc[i, 'SC'])
    end = time.time()
    print((end-start)/60)

    pccs = np.corrcoef(np.array(exp_growth),
                       np.array(growth))
    # pairp = stats.ttest_rel(np.array(exp.loc[:100, 'SC'].values.tolist()),
    #                         np.array(growth)
    #                         , axis=0, nan_policy='propagate')
    print('PCC:' + str(pccs[0, 1]))


    font1 = {'family': 'Arial',
             'weight': 'normal',
             'size': 23,
             }
    font2 = {'family': 'Arial',
             'weight': 'normal',
             'size': 18,
             }
    font3 = {'family': 'Arial',
             'weight': 'normal',
             'size': 12,
             }
    plt.scatter(x = np.array(exp_growth),
                y = np.array(growth),
                alpha=0.3)
    plt.ylabel('simulation',
               fontdict=font2)
    plt.xlabel('Growth rate',
               fontdict=font2)
    plt.xticks(fontproperties='Arial',
               size=16)
    plt.yticks(fontproperties='Arial',
               size=16)
    plt.tight_layout()
    # plt.savefig('../deletion/output/prot_chengfa.jpg')
    plt.show()
