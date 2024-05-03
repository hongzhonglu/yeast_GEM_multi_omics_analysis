###########################################################################
# run this code and get regulation analysis histogram
###########################################################################
import cobra
import pandas as pd
from N_lim import *
from N_lim.code.getflux_pro import *
from N_lim.code.load_csmodel import load_csmodel
import numpy as np
import scipy.stats as st
from N_lim.code.get_coefficient import *
from N_lim.code.pairing import paring
from N_lim.code.Histogram import histogram
import os

os.chdir('..')


def get_flux_prot(ID, model_fluxdict, protdict):
    cc = 0
    reaction_slope = pd.DataFrame()
    for rea in ID:
        flux_prot = paring(model_fluxdict, protdict, rea)
        slope = get_coefficient(flux_prot)
        if slope != 'nan':
            reaction_slope.loc[rea, 'p'] = slope
        if np.sum(flux_prot.loc[:, 'flux']) > 0.000001:
            cc += 1
    print(cc)
    c = 0
    for r in reaction_slope.loc[:, 'p']:
        if r != 'nan':
            c += 1
    print(c)
    return reaction_slope, cc, c


def get_flux_prot_z(ID, model_fluxdict, protdict):
    cc = 0
    reaction_slope = pd.DataFrame()
    for rea in ID:
        flux_prot = paring(model_fluxdict, protdict, rea)
        slope = get_coefficient(flux_prot)
        if slope != 'nan':
            reaction_slope.loc[rea, 'p'] = slope
        if np.sum(flux_prot.loc[:, 'flux']) > 0.000001:
            cc += 1
    print(cc)
    c = 0
    for r in reaction_slope.loc[:, 'p']:
        if r != 'nan':
            c += 1
    print(c)
    return reaction_slope, cc, c


def get_id(sub):
    y9 = pd.read_excel("../N_lim/data/yeast-GEM.xlsx")
    y9.index = y9.loc[:, 'ID']
    rxn = [r for r in y9.index if y9.loc[r, 'subsystem'] in sub]
    return rxn


def z_score(data):
    std_data = {}
    for k, v in data.items():
        new = v.copy()
        val = [(i-np.mean(v.iloc[:, 1]))/np.std(v.iloc[:, 1]) for i in v.iloc[:, 1]]
        new.iloc[:, 1] = val
        std_data[k] = new
    return std_data


[modelCN115, modelCN50, modelCN30, modelphe_Nlim_01,
modelile_Nlim_01, modelNH4_Nlim_035, modelNH4_Nlim_030,
modelNH4_Nlim_018, modelNH4_Nlim_013, modelNH4_Nlim_01,
modelNH4_Nlim_005, modelgln_Nlim_01] = load_csmodel()
modelname = ['modelCN115', 'modelCN50', 'modelCN30', 'modelphe_Nlim_01',
             'modelile_Nlim_01', 'modelNH4_Nlim_035', 'modelNH4_Nlim_030',
             'modelNH4_Nlim_018', 'modelNH4_Nlim_013', 'modelNH4_Nlim_01',
             'modelNH4_Nlim_005', 'modelgln_Nlim_01']
modellist = [modelCN115, modelCN50, modelCN30, modelphe_Nlim_01,
             modelile_Nlim_01, modelNH4_Nlim_035, modelNH4_Nlim_030,
             modelNH4_Nlim_018, modelNH4_Nlim_013, modelNH4_Nlim_01,
             modelNH4_Nlim_005, modelgln_Nlim_01]
conditionname = ['CN115', 'CN50', 'CN30', 'Phe', 'Ile',
                 'N035', 'N030', 'N018', 'N013',
                 'N010', 'N005', 'Gln']

# choose pFBA or sampling
allg = 'pFBA'

# get flux and protein data
model_fluxdict = {}
condition_prodict = {}
print('compute flux')
for m, mo in zip(modelname, modellist):
    print(m)
    tempflux = getflux(mo, allg, m)
    model_fluxdict[m] = tempflux
print('compute protein')
for c in conditionname:
    tempprot = getpro(c)
    condition_prodict[c] = tempprot
for r in model_fluxdict.values():
    r.index = r.loc[:, 'reactionID']
for co in condition_prodict.values():
    co.index = co.loc[:, 'gene']

# get commmon protein of two dataset
N_protgene_list = condition_prodict['Phe'].index.values.tolist()
CN_protgene_list = condition_prodict['CN30'].index.values.tolist()
commongene = [val for val in N_protgene_list if val in CN_protgene_list]
protdict = {}
for con, ke in zip(condition_prodict.values(), condition_prodict.keys()):
    prottemp = {}
    for v in commongene:
        prottemp[v] = con.loc[v, 'prot']
    temp = pd.DataFrame(index=prottemp.keys(),
                        columns=['gene', 'prot'])
    temp.loc[:, 'gene'] = commongene
    temp.loc[:, 'prot'] = prottemp.values()
    protdict[ke] = temp


# z-score
# protdict = z_score(protdict)
# model_fluxdict = z_score(model_fluxdict)

# compute slope
ID = model_fluxdict['modelCN115'].index.values.tolist()
# all_reaction_slope, all_cc, all_c = get_flux_prot(ID, model_fluxdict, protdict)
all_reaction_slope, all_cc, all_c = get_flux_prot_z(ID, model_fluxdict, protdict)

morethan05 = [r for r in all_reaction_slope.index if all_reaction_slope.loc[r, 'p'] > 0.5]
m05 = pd.DataFrame(morethan05)
m05.to_csv(r'./output/m05.csv')

carbon_rxn = get_id(['Citrate cycle (TCA cycle)',
                      'Glycolysis / gluconeogenesis',
                      'Pentose phosphate pathway',
                      'Fructose and mannose metabolism',
                      'Galactose metabolism',
                      'Ascorbate and aldarate metabolism',
                      'Starch and sucrose metabolism',
                      'Amino sugar and nucleotide sugar metabolism',
                      'Pyruvate metabolism',
                      'Glyoxylate and dicarboxylate metabolism',
                      'Propanoate metabolism',
                      'Butanoate metabolism',
                      'C5-Branched dibasic acid metabolism',
                      'Inositol phosphate metabolism'])
carbon_reaction_slope, carbon_cc, carbon_c = get_flux_prot(carbon_rxn, model_fluxdict, protdict)
amino_rxn = get_id(['Alanine, aspartate and glutamate metabolism',
                    'Glycine, serine and threonine metabolism',
                    'Cysteine and methionine metabolism',
                    'Valine, leucine and isoleucine metabolism',
                    'Lysine metabolism',
                    'Arginine biosynthesis',
                    'Arginine and proline metabolism',
                    'Histidine metabolism',
                    'Tyrosine metabolism',
                    'Phenylalanine metabolism',
                    'Tryptophan metabolism',
                    'Phenylalanine, tyrosine and tryptophan biosynthesis'])
amino_reaction_slope, amino_cc, amino_c = get_flux_prot(amino_rxn, model_fluxdict, protdict)
lipid_rxn = get_id(['Glycerophospholipid metabolism',
                    'Biosynthesis of unsaturated fatty acids',
                    'Glycerolipid metabolism',
                    'Steroid biosynthesis',
                    'Fatty acid biosynthesis',
                    'Sphingolipid metabolism',
                    'Arachidonic acid metabolism'])
lipid_reaction_slope, lipid_cc, lipid_c = get_flux_prot(lipid_rxn, model_fluxdict, protdict)


reaction_slope = [all_reaction_slope, carbon_reaction_slope,
                  amino_reaction_slope]
histogram(reaction_slope)
