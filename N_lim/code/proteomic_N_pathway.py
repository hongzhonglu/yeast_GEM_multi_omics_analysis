# 找到yeastGEM中以NH4为底物的反应，重点分析其基因， 随后次要分析kegg的subsystem

import pandas as pd
import numpy as np
import cobra


class Gene_sub(cobra.Gene):
    def __init__(self, name, sub):
        super().__init__()
        self.subsystem = sub
        self.name = name


def find_ammonium(rxn):
    if rxn.reversibility == False:
        name = [m.name for m in rxn.reactants if m.name == 'ammonium']
    else:
        name = [m.name for m in rxn.reactants + rxn.products if m.name == 'ammonium']
    if len(name) == 1:
        cont_N = rxn
    else:
        cont_N = False
    return cont_N


def get_gene_name(gene):
    import re
    g = re.compile(r'Y\w{6}')
    out = g.findall(gene)
    return out


def get_gene_sub(n_gene):
    sub_dic = dict()
    for gene, sub in zip(n_gene.loc[:, 'gpr'], n_gene.loc[:, 'sub']):
        gene_list = get_gene_name(gene)
        for g in gene_list:
            if g in sub_dic.keys():
                sub_dic[g] = sub_dic[g] + '+' + sub
            else:
                sub_dic[g] = sub
    gene_sub = {}
    for gene, sub in zip(sub_dic.keys(), sub_dic.values()):
        gene_sub[gene] = list(set(sub.split('+')))
    return gene_sub


def std(data, dirc):
    new = data.copy()
    if dirc == 1:
        for i in data.index:
            rawdata = data.loc[i, ['NH4', 'Gln', 'Phe', 'Ile']]
            # newdata = [(x-np.mean(rawdata))/np.std(rawdata) for x in rawdata]
            # newdata = [(x - np.min(rawdata)) / (np.max(rawdata)-np.min(rawdata))for x in rawdata]
            newdata = [x / (np.sum(rawdata)) for x in rawdata]
            new.loc[i, ['NH4', 'Gln', 'Phe', 'Ile']] = newdata
    if dirc == 0:
        for i in data.columns:
            rawdata = data.loc[:, i]
            # newdata = [(x-np.mean(rawdata))/np.std(rawdata) for x in rawdata]
            # newdata = [(x - np.min(rawdata)) / (np.max(rawdata)-np.min(rawdata))for x in rawdata]
            newdata = [x / (np.sum(rawdata)) for x in rawdata]
            new.loc[:, i] = newdata
    return new


def log2fc(data):
    new = data.copy()
    for i in data.index:
        rawdata = list(data.loc[i, ['NH4', 'Gln', 'Phe', 'Ile']])
        newdata = [x / rawdata[0] for x in rawdata]
        new.loc[i, ['NH4', 'Gln', 'Phe', 'Ile']] = np.log2(newdata)
    return new


def get_sub_gene(sub, sgd, prot_N, gene_sub):
    ks = list(prot_N.index)
    vs = [gene_sub[sgd.loc[k, 'Systematic_name']] for k in ks]
    genes = [str(k) for k, v in zip(ks, vs) if sub in v]
    genes = list(set(genes))
    return genes


def comp_sub_resource(sub_gene, prot_N):
    result = pd.DataFrame(index=sub_gene.keys(), columns=['NH4', 'Gln', 'Phe', 'Ile'])
    for k, v in zip(sub_gene.keys(), sub_gene.values()):
        nh4 = []
        gln = []
        phe = []
        ile = []
        for vv in v:
            nh4.append(prot_N.loc[vv, 'NH4'])
            gln.append(prot_N.loc[vv, 'Gln'])
            phe.append(prot_N.loc[vv, 'Phe'])
            ile.append(prot_N.loc[vv, 'Ile'])
        nh4_conc = np.sum(nh4)
        gln_conc = np.sum(gln)
        phe_conc = np.sum(phe)
        ile_conc = np.sum(ile)
        result.loc[k, :] = [nh4_conc, gln_conc, phe_conc, ile_conc]
    return result


model = cobra.io.read_sbml_model(r'C:\Users\yuhuzhouye\Desktop\yeast9_w\yesat9_temp\yeast-GEM.xml')
n_gene = pd.read_csv('../output/N_rxns.csv')
prot = pd.read_csv('../data/proteomic.csv')
sgd = pd.read_table('../data/SGDgeneNames.tsv')
all_sub = target_sub = [
    'Purine metabolism',
    'Pyrimidine metabolism',
    'Alanine, aspartate and glutamate metabolism',
    'Glycine, serine and threonine metabolism',
    'Cysteine and methionine metabolism',
    'Valine, leucine and isoleucine degradation',
    'Valine, leucine and isoleucine biosynthesis',
    'Lysine biosynthesis',
    'Lysine degradation',
    'Arginine biosynthesis',
    'Arginine and proline metabolism',
    'Histidine metabolism',
    'Tyrosine metabolism',
    'Phenylalanine metabolism',
    'Tryptophan metabolism',
    'Phenylalanine, tyrosine and tryptophan biosynthesis',
]

# find all gene related to nitrogen usage
# rxn_ammonium = [r for r in model.reactions if find_ammonium(r) != False]
# gene_ammonium = [r.gene_reaction_rule for r in rxn_ammonium]
all_gene = ''
for g in list(n_gene.loc[:, 'gpr']):
    all_gene = all_gene + g
all_gene = list(set(get_gene_name(all_gene)))
gene_sub = get_gene_sub(n_gene)

# get target prot conc
prot_N = pd.DataFrame(columns=prot.columns)
sgd.index = sgd.loc[:, 'Standard_name']
prot.index = prot.loc[:, 'Gene']
for ge in prot.index:
    try:
        if sgd.loc[ge, 'Systematic_name'] in gene_sub.keys():
            prot_N.loc[ge, :] = prot.loc[ge, :]
    except KeyError:
        pass


# get sub gene
sub_gene = {k: get_sub_gene(k, sgd, prot_N, gene_sub) for k in all_sub}
# missing genes related with isoleucine are from kegg
# CHA1, ALD3, LEU9, POT1 are not in prot
sub_gene['Valine, leucine and isoleucine degradation'] = ['BAT1', 'BAT2', 'LPD1',
                                                          'EHD3', 'ERG13',
                                                          'ERG10', 'UGA1', 'ALD2',
                                                          'ALD4', 'ALD5',
                                                          'ALD6', 'HFD1']
sub_gene['Valine, leucine and isoleucine biosynthesis'] = ['LEU2', 'ILV1',
                                                           'ILV2', 'ILV6', 'ILV5',
                                                           'ILV3', 'BAT1', 'BAT2',
                                                           'LEU1', 'LEU4',
                                                           ]
sub_gene['Lysine degradation'] = ['LYS1', 'LYS9', 'KGD2', 'LPD1', 'UGA2',
                                  'ERG10', 'SET1', 'SET2',
                                ]
sub_gene['Lysine biosynthesis'] = ['LYS21', 'LYS20', 'LYS4', 'LYS12', 'ARO8',
                                   'LYS2', 'LYS9', 'LYS1',
                                ]
gene_from_kegg = list(set(sub_gene['Valine, leucine and isoleucine degradation'] +
                          sub_gene['Valine, leucine and isoleucine biosynthesis'] +
                          sub_gene['Lysine biosynthesis'] +
                          sub_gene['Lysine degradation']))
for g in gene_from_kegg:
    try:
        prot_N.loc[g, :] = prot.loc[g, :]
    except KeyError:
        print(g)


# compute fc
# std_gene = std(prot_N)
# log2fc_gene = log2fc(std_gene)
# std_gene.drop(['Accession', 'Gene'], axis = 1, inplace=True)
# std_gene.to_excel('../output/std_gene.xlsx')


final_N_resource = comp_sub_resource(sub_gene, prot_N)
final_N_resource.to_excel('../output/raw_pathw.xlsx')
# std_pathw = std(final_N_resource, dirc=0)
std_pathwh = std(final_N_resource, dirc=1)
log2fc_pathw = log2fc(std_pathwh)
std_pathwh.to_excel('../output/std_pathwh.xlsx')
std_pathwl = std(final_N_resource, dirc=0)
log2fc_pathw = log2fc(std_pathwl)
std_pathwl.to_excel('../output/std_pathwl.xlsx')

