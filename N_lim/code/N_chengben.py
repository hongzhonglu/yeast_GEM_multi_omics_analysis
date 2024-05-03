# Define a parameter, preference score, to describe the preference of  yeast using a nitrogen source.
#
# How to compute:
# First, set nitrogen source exchange rate as objective, and maximize it to get minimum exchange rate(NME).
# Set glucose exchange rate as objective.
# Change nitrogen source exchange rate in a small scale( 1~1.5 times NME) and compute maximal glucose exchange rate.
# Get absolute value of slope.
# get "Preference Score.jpg"

import cobra
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.pyplot import MultipleLocator

os.chdir('..')

def N_chengben(model, rxn):
    mini_b = model.medium[rxn]
    intake = np.linspace(mini_b, 1.5*mini_b, 5)
    atp = []
    for i in intake:
        tempmedium = model.medium
        tempmedium[rxn] = i
        tempmedium['r_1714'] = 1000
        model.medium = tempmedium
        model.reactions.get_by_id('r_1714').bounds = (-1000, 0)
        model.reactions.get_by_id('r_2111').bounds = (0.1, 0.1)
        model.objective = model.reactions.get_by_id('r_1714')
        model.objective_direction = 'max'
        solu = model.optimize()
        atp.append(solu.objective_value)
    slope, intercept = np.polyfit(intake, atp, 1)
    return slope


# def themrConstrain(model):
#     # themrBack = pd.read_csv(r'../../thermodynamic analysis/output/thermobackw.csv', index_col=False)
#     themrBack = pd.read_csv(r"D:\model_research\yeast_GEM_multi_omics_analysis\thermodynamic analysis\output\thermobackw.csv", index_col=False)
#     themrBack = themrBack.columns.tolist()
#     themrforward = pd.read_csv(r"D:\model_research\yeast_GEM_multi_omics_analysis\thermodynamic analysis\output\thermoforw.csv", index_col=False)
#     themrforward = themrforward.columns.tolist()
#     themrbi = pd.read_csv(r"D:\model_research\yeast_GEM_multi_omics_analysis\thermodynamic analysis\output\thermobi.csv", index_col=False)
#     themrbi = themrbi.columns.tolist()
#     for r in model.reactions:
#         for b in themrBack:
#             if r.id == b:
#                 model.reactions.get_by_id(b).bounds = (-1000, 0)
#         for f in themrforward:
#             if r.id == f:
#                 model.reactions.get_by_id(f).bounds = (0, 1000)
#         for bi in themrbi:
#             if r.id == bi:
#                 model.reactions.get_by_id(bi).bounds = (-1000, 1000)
#     # fix NGAM
#     model.reactions.get_by_id('r_4046').bounds = (0.7, 0.7)
#     # fix glucose intake
#     model.reactions.get_by_id('r_1714').bounds = (-1, 1000)
#
#     return model


def StandardScaler(data):
    StandardData = [i/np.sum(data) for i in data]
    return StandardData


path = [
    '../N_lim/output/Nlim_model/modelgln_Nlim_01.xml',
    '../N_lim/output/Nlim_model/modelNH4_Nlim_01.xml',
    '../N_lim/output/Nlim_model/modelphe_Nlim_01.xml',
    '../N_lim/output/Nlim_model/modelile_Nlim_01.xml',
]
rx = [
    'r_1891',
    'r_1654',
    'r_1903',
    'r_1897',
]
abssl = []
for m, r in zip(path, rx):
    model = cobra.io.read_sbml_model(m)
    # medium = model.medium
    # #model = themrConstrain(model)
    # model.medium = medium
    slope = N_chengben(model, r)
    abssl.append(abs(slope))
abssl = StandardScaler(abssl)

font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 28,
         }
font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 24,
         }
font3 = {'family': 'Arial',
         'weight': 'normal',
         'size': 20,
         }

plt.figure(figsize=(4, 6))
plt.bar(['Glutamine', 'Ammonium', 'Phenylalanine', 'Isoleucine'],
        abssl,
        facecolor='#C69C6D',
        edgecolor='#000000')
plt.ylabel('Scaled preference score', fontdict=font2)

y_major_locator = MultipleLocator(1)
ax = plt.gca()
#ax.yaxis.set_major_locator(y_major_locator)
plt.yticks(fontproperties='Arial',
           size=18)
plt.xticks(fontproperties='Arial',
           size=18,
           rotation=30)
plt.ylim(ymax=0.38)
plt.ylim(ymin=0.15)

plt.tight_layout()

plt.savefig('../N_lim/output/Preference score.jpg',
            dpi=600)
plt.show()

