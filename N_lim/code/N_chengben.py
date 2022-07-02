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
        model.objective = model.reactions.get_by_id('r_1714')
        model.objective_direction = 'min'
        solu = model.optimize()
        atp.append(solu.objective_value)
    slope, intercept = np.polyfit(intake, atp, 1)
    return slope

path = [
    '../N_lim/output/Nlim_model/modelgln_Nlim_01.xml',
    '../N_lim/output/Nlim_model/modelNH4_Nlim_01.xml',
    '../N_lim/output/Nlim_model/modelile_Nlim_01.xml',
    '../N_lim/output/Nlim_model/modelphe_Nlim_01.xml',
]
rx = [
    'r_1891',
    'r_1654',
    'r_1897',
    'r_1903',
]
abssl = []
for m, r in zip(path, rx):
    model = cobra.io.read_sbml_model(m)
    slope = N_chengben(model, r)
    abssl.append(abs(slope))

font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 23,
         }
font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 16,
         }
font3 = {'family': 'Arial',
         'weight': 'normal',
         'size': 12,
         }

plt.subplots(figsize=(8, 6))
plt.bar(['glutamine', 'ammonium', 'isoleucine', 'phenylalanine'],
        abssl)
plt.xlabel('Nitrogen Source', fontdict=font2)
plt.ylabel('Preference Score', fontdict=font2)
plt.xticks(fontproperties='Arial',
           size=12)
plt.yticks(fontproperties='Arial',
           size=12)
plt.savefig('../N_lim/output/Preference Score.jpg',
            dpi=600)
plt.show()
