##########################################
# get "Nsourceflux.xlsx'"


import cobra
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

os.chdir('..')

path = [
    '../N_lim/output/Nlim_model/modelgln_Nlim_01.xml',
    '../N_lim/output/Nlim_model/output/modelNH4_Nlim_01.xml',
    '../N_lim/output/Nlim_model/output/modelile_Nlim_01.xml',
    '../N_lim/output/Nlim_model/output/modelphe_Nlim_01.xml',
]
rx = [
    'r_1891',
    'r_1654',
    'r_1897',
    'r_1903',
]
model = cobra.io.read_sbml_model(
    '../N_lim/output/Nlim_model/modelgln_Nlim_01.xml')
flux = pd.DataFrame(index=model.optimize().fluxes.index, columns=rx)
for m, r in zip(path, rx):
    model = cobra.io.read_sbml_model(m)
    tempmedium = model.medium
    tempmedium[r] = 1000
    model.medium = tempmedium
    model.reactions.get_by_id(r).bounds = (-1000, 0)
    model.objective_direction = 'max'
    so = model.optimize()
    flux.loc[:, r] = so.fluxes[:]

flux2 = flux
flux2.columns = ['Glu', 'NH4', 'Ile', 'Phe']
flux2.to_excel('../N_lim/output/Nlim_model/output/Nsourceflux.xlsx')
def tonglu(flux, metabolism):
    subfile = pd.read_table(
        '../N_lim/data/uniqueSubsystems.tsv')
    rxn = []
    for i in subfile.index.values.tolist():
        if subfile.loc[i, 'subsystem_unique_@2021_12'] == metabolism:
            rxn.append(subfile.loc[i, 'ID'])
    tltemp = []
    tl = []
    for n in flux.columns.values.tolist():
        tempflux = 0
        co = 0
        for mr in rxn:
            if flux.loc[mr, n] > 0.000001:
                tempflux += abs(flux.loc[mr, n])
                co += 1
        meanflux = tempflux / co
        tltemp.append(meanflux)

    if np.max(tltemp) - np.min(tltemp) != 0:
        tl = [((st - np.min(tltemp)) / (np.max(tltemp) - np.min(tltemp))) for st in tltemp]
    else:
        print(metabolism + ': zero')
    return tl

# Carbohydrate metabolism
c_gly = tonglu(flux, 'Glycolysis / gluconeogenesis')
c_TCA = tonglu(flux, 'Citrate cycle (TCA cycle)')
c_ppp = tonglu(flux, 'Pentose phosphate pathway')
# c_pgi = tonglu(flux, 'Pentose and glucuronate interconversions')
# c_fmm = tonglu(flux, 'Fructose and mannose metabolism')
c_gm = tonglu(flux, 'Galactose metabolism')
# c_aam = tonglu(flux, 'Ascorbate and aldarate metabolism')
c_pm = tonglu(flux, 'Pyruvate metabolism')

# Energy metabolism
e_op = tonglu(flux, 'Oxidative phosphorylation')
# e_mm = tonglu(flux, 'Methane metabolism')
e_nm = tonglu(flux, 'Nitrogen metabolism')
e_sm = tonglu(flux, 'Sulfur metabolism')

# Amino acid metabolism
a_sm = tonglu(flux, 'Alanine, aspartate and glutamate metabolism')
a_gstm = tonglu(flux, 'Glycine, serine and threonine metabolism')
a_cmm = tonglu(flux, 'Cysteine and methionine metabolism')
a_vli = tonglu(flux, 'Valine, leucine and isoleucine metabolism')
a_lm = tonglu(flux, 'Lysine metabolism')
a_ab = tonglu(flux, 'Arginine biosynthesis')
a_apm = tonglu(flux, 'Arginine and proline metabolism')
a_hm = tonglu(flux, 'Histidine metabolism')
# a_tym = tonglu(flux, 'Tyrosine metabolism')
# a_pm = tonglu(flux, 'Phenylalanine metabolism')
a_tm = tonglu(flux, 'Tryptophan metabolism')
a_ptt = tonglu(flux, 'Phenylalanine, tyrosine and tryptophan biosynthesis')

# Nucleotide metabolism
n_pur = tonglu(flux, 'Purine metabolism')
n_pyr = tonglu(flux, 'Pyrimidine metabolism')

pathway = np.vstack((c_gly, c_TCA, c_ppp, c_gm, c_pm,
                     e_op, e_nm, e_sm,
                     a_sm, a_gstm, a_cmm, a_vli, a_lm, a_ab, a_apm, a_hm, a_tm, a_ptt,
                     n_pur, n_pyr))
pathway = pd.DataFrame(pathway,
                       index=['Glycolysis / gluconeogenesis', 'Citrate cycle (TCA cycle)', 'Pentose phosphate pathway',
                              'Galactose metabolism', 'Pyruvate metabolism',
                              'Oxidative phosphorylation', 'Nitrogen metabolism', 'Sulfur metabolism',
                              'Alanine, aspartate and glutamate metabolism', 'Glycine, serine and threonine metabolism',
                              'Cysteine and methionine metabolism', 'Valine, leucine and isoleucine metabolism',
                              'Lysine metabolism', 'Arginine biosynthesis', 'Arginine and proline metabolism',
                              'Histidine metabolism', 'Tryptophan metabolism',
                              'Phenylalanine, tyrosine and tryptophan biosynthesis', 'Purine metabolism',
                              'Pyrimidine metabolism'],
                       columns=['glutamine', 'ammonium', 'isoleucine', 'phenylalanine'])
f, ax = plt.subplots(figsize=(10, 8))
h = sns.heatmap(pathway,
                cmap="YlGnBu")
# plt.xticks(fontproperties='Arial',
#            size=12,
#            )
# plt.yticks(fontproperties='Arial',
#            size=12,
#            )
# 把横纵轴展示全
plt.tight_layout()
plt.savefig('../N_lim/output/pathAna.jpg',
            dpi=600)
plt.show()

###################################
