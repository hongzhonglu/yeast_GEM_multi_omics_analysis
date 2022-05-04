###########################################################################
# run this code and get histogram
###########################################################################
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
           'all\n'
           'central carbon metabolism\n'
           'amino')
ID = [d for d in id.split(' ')]

if id == 'all':
    ID = model_fluxdict['modelCN115'].index.values.tolist()
if id =='central carbon metabolism':
    ID = ['r_0091', 'r_0091', 'r_0163', 'r_0165', 'r_0173',
          'r_0174', 'r_0280', 'r_0300', 'r_0301', 'r_0302',
          'r_0303', 'r_0322', 'r_0332', 'r_0356', 'r_0366',
          'r_0449', 'r_0450', 'r_0451', 'r_0452', 'r_0466',
          'r_0467', 'r_0486', 'r_0533', 'r_0534', 'r_0535',
          'r_0658', 'r_0659', 'r_0661', 'r_0713', 'r_0714',
          'r_0715', 'r_0831', 'r_0832', 'r_0884', 'r_0886',
          'r_0887', 'r_0888', 'r_0889', 'r_0892', 'r_0893',
          'r_0907', 'r_0916', 'r_0959', 'r_0960', 'r_0961',
          'r_0962', 'r_0969', 'r_0982', 'r_0984', 'r_0990',
          'r_1021', 'r_1022', 'r_1048', 'r_1049', 'r_1050',
          'r_1054', 'r_2115', 'r_2116', 'r_2131', 'r_2305',
          'r_4232', 'r_4283', 'r_4284', 'r_4285', 'r_4286',
          'r_4287', 'r_4288', 'r_4566'
         ]
if id == 'amino':
    ID = ['r_0007', 'r_0012', 'r_0016', 'r_0018', 'r_0023', 'r_0024', 'r_0025', 'r_0027',
          'r_0029', 'r_0030', 'r_0045', 'r_0058', 'r_0060', 'r_0061', 'r_0068', 'r_0096',
          'r_0097', 'r_0145', 'r_0156', 'r_0176', 'r_0177', 'r_0178', 'r_0186', 'r_0187',
          'r_0199', 'r_0200', 'r_0201', 'r_0211', 'r_0215', 'r_0216', 'r_0217', 'r_0218',
          'r_0219', 'r_0225', 'r_0309', 'r_0352', 'r_0353', 'r_0443', 'r_0468', 'r_0469',
          'r_0470', 'r_0471', 'r_0472', 'r_0473', 'r_0475', 'r_0476', 'r_0477', 'r_0500',
          'r_0501', 'r_0502', 'r_0503', 'r_0504', 'r_0505', 'r_0506', 'r_0507', 'r_0508',
          'r_0509', 'r_0536', 'r_0537', 'r_0538', 'r_0541', 'r_0542', 'r_0545', 'r_0546',
          'r_0547', 'r_0548', 'r_0563', 'r_0564', 'r_0567', 'r_0663', 'r_0664', 'r_0669',
          'r_0670', 'r_0671', 'r_0672', 'r_0673', 'r_0674', 'r_0675', 'r_0676', 'r_0678',
          'r_0679', 'r_0680', 'r_0681', 'r_0682', 'r_0683', 'r_0687', 'r_0689', 'r_0690',
          'r_0692', 'r_0693', 'r_0694', 'r_0699', 'r_0700', 'r_0762', 'r_0763', 'r_0817',
          'r_0819', 'r_0891', 'r_0909', 'r_0910', 'r_0917', 'r_0918', 'r_0920', 'r_0929',
          'r_0936', 'r_0937', 'r_0940', 'r_0957', 'r_0988', 'r_0989', 'r_1001', 'r_1002',
          'r_1023', 'r_1040', 'r_1041', 'r_1063', 'r_1065', 'r_1087', 'r_1088', 'r_2112',
          'r_2113', 'r_2114', 'r_4042', 'r_1739', 'r_1871', 'r_1872', 'r_1887', 'r_2050',
          'r_4192', 'r_4193', 'r_4223', 'r_4226', 'r_4228', 'r_4248', 'r_4254', 'r_4255',
          'r_4262', 'r_4273', 'r_4292', 'r_4293', 'r_4294', 'r_4298', 'r_4299', 'r_4301',
          'r_4302', 'r_4312', 'r_4336', 'r_4363', 'r_4366', 'r_4407', 'r_4485', 'r_4486',
          'r_4487', 'r_4488', 'r_4570', 'r_4576', 'r_4579', 'r_4581',
    ]
# compute slope
reaction_slope = {}
for rea in ID:
    flux_prot = paring(model_fluxdict, protdict, rea)
    slope = get_coefficient(flux_prot)
    reaction_slope[rea] = slope

name = input('what is your name? i mean the figure\n'
             'for example\n'
             'my_figure_name_for_sample')
# get graph
histogram(reaction_slope, savename=name)

