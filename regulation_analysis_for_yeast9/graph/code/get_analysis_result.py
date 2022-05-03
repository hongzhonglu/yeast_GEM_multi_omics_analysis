import cobra
import pandas as pd
from regulation_analysis_for_yeast9 import *
from regulation_analysis_for_yeast9.get_fluxandproteoimcs.code.getflux_pro import *
from regulation_analysis_for_yeast9.get_fluxandproteoimcs.code.load_csmodel import load_csmodel
import numpy as np
import scipy.stats as st
from regulation_analysis_for_yeast9.get_regulation_analysis_coefficient.code.get_coefficient import *
from regulation_analysis_for_yeast9.get_regulation_analysis_coefficient.code.pairing import paring
from regulation_analysis_for_yeast9.graph.code.Histogram import histogram

# load condition specific model modified in matlab

[modelphe_Nlim_01, modelile_Nlim_01, modelNH4_Nlim_035,
 modelNH4_Nlim_030, modelNH4_Nlim_018, modelNH4_Nlim_013,
 modelCN115, modelCN50, modelCN30, modelNH4_Nlim_01,
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
allg = input('pFBA or sample (⑅˃◡˂⑅)')

# get flux and protein data
model_fluxdict = {}
condition_prodict = {}
for m, mo in zip(modelname, modellist):
    tempflux = getflux(mo, allg)
    model_fluxdict[m] = tempflux
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
    temp.loc[:,'gene'] = commongene
    temp.loc[:,'prot'] = prottemp.values()
    protdict[ke] = temp

# choose reactions
id = input('please us your reaction list (⑅˃◡˂⑅)\nfor example\n'
           'r_0001\n'
           'r_0001 r_0002 r_0003\n'
           'all')
ID = [d for d in id.split(' ')]

if id == 'all':
    ID = model_fluxdict['modelCN115'].index.values.tolist()

# compute slope
reaction_slope = {}
for rea in ID:
    flux_prot = paring(model_fluxdict, protdict, rea)
    slope = get_coefficient(flux_prot)
    reaction_slope[rea] = slope
# get graph
histogram(reaction_slope)

