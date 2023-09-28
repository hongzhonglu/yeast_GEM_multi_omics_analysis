##########################################
# get "ile.jpg","phe.jpg","phe/ile_zgene.xlsx"
# first run "FluxAnalysis.py" and get "Nsourceflux.xlsx"
# Select the reaction whose flux and corresponding protein expression level of
# Ile or Phe that are larger than their counterpart in NH4.
# Isoenzyme: sum of expression level
# Complex: minimum of all subunits

import pandas as pd
import cobra
import numpy as np
import matplotlib.pyplot as plt
import os




def onetoone(rxnid, gr, Nsource, flux, pro):
    zid = 'nan'
    fid = 'nan'
    genename = 'nan'
    try:
        genename = SGD.loc[gr, 'Standard_name']
        fluxdiff = flux.loc[rxnid, Nsource] - flux.loc[rxnid, 'NH4']
        prodiff = pro.loc[genename, Nsource] - pro.loc[genename, 'NH4']
        pro_Nsource = pro.loc[genename, Nsource]
        pro_NH4 = pro.loc[genename, 'NH4']
        if fluxdiff > 0 and prodiff > 0:
            zid = rxnid
        if fluxdiff < 0 and prodiff < 0:
            fid = rxnid
    except KeyError:
        pro_Nsource = 'nan'
        pro_NH4 = 'nan'
    return zid, fid, pro_Nsource, pro_NH4, genename

# sum
def isoenzyme(o_s, rxnid, Nsource, flux, pro, SGD):
    zid = 'nan'
    fid = 'nan'
    genename = []
    fluxdiff = flux.loc[rxnid, Nsource] - flux.loc[rxnid, 'NH4']
    try:
        # compare flux and proteomic data
        pro_N_l = [pro.loc[SGD.loc[sn, 'Standard_name'], Nsource] for sn in o_s]
        pro_NH4_l = [pro.loc[SGD.loc[sn, 'Standard_name'], 'NH4'] for sn in o_s]
        pro_Nsource = np.sum(pro_N_l)
        pro_NH4 = np.sum(pro_NH4_l)
        prodiff = pro_Nsource - pro_NH4
        if fluxdiff > 0 and prodiff > 0:
            # get rxn id
            zid = rxnid
            # get gene id
            for v1, v2 in zip(pro_N_l, pro_NH4_l):
                if v1 > v2:
                    genename.append(o_s[pro_N_l.index(v1)])
        if fluxdiff < 0 and prodiff < 0:
            fid = rxnid

    except KeyError:
        pass
    return zid, fid, genename

# min
def complex(gr, rxnid, Nsource, flux, pro, SGD):
    gr_s = gr.split(' and ')
    zid = 'nan'
    fid = 'nan'
    fluxdiff = flux.loc[rxnid, Nsource] - flux.loc[rxnid, 'NH4']
    genename = []
    try:
        pro_N_l = [pro.loc[SGD.loc[sn, 'Standard_name'], Nsource] for sn in gr_s]
        pro_NH4_l = [pro.loc[SGD.loc[sn, 'Standard_name'], 'NH4'] for sn in gr_s]
        pro_Nsource = np.min(pro_N_l)
        pro_NH4 = np.min(pro_NH4_l)
        prodiff = pro_Nsource - pro_NH4
        if fluxdiff > 0 and prodiff > 0:
            zid = rxnid
            for v1, v2 in zip(pro_N_l, pro_NH4_l):
                if v1 > v2:
                    genename.append(gr_s[pro_N_l.index(v1)])
        if fluxdiff < 0 and prodiff < 0:
            fid = rxnid
    except KeyError:
        pro_Nsource = 'nan'
        pro_NH4 = 'nan'
    return zid, fid, pro_Nsource, pro_NH4, genename


def comiso(o_s, rxnid, Nsource, flux, pro, SGD):
    zid = 'nan'
    fid = 'nan'
    genename = []
    tempcom = []
    tempiso = []
    fluxdiff = flux.loc[rxnid, Nsource] - flux.loc[rxnid, 'NH4']
    for i in o_s:
        t = i.split(' and ')
        if len(t) == 1:
            tempiso.append(i)
        if len(t) != 1:
            t1 = i.strip('( ')
            t2 = t1.strip(' )')
            tempcom.append(t2)
    xiaoNsource = []
    xiaoNH4 = []
    for co in tempcom:
        zid, fid, pro_Nsource, pro_NH4, genenamecom = complex(co, rxnid, Nsource, flux, pro, SGD)
        for gi in genenamecom:
            genename.append(gi)
        if pro_Nsource != 'nan' and pro_NH4 != 'nan':
            xiaoNsource.append(pro_Nsource)
            xiaoNH4.append(pro_NH4)
    for iso in tempiso:
        zid, fid, pro_Nsource, pro_NH4, genenameone = onetoone(rxnid, iso, Nsource, flux, pro)
        genename.append(genenameone)
        if pro_Nsource != 'nan' and pro_NH4 != 'nan':
            xiaoNsource.append(pro_Nsource)
            xiaoNH4.append(pro_NH4)
    prodiff = np.sum(xiaoNsource) - np.sum(xiaoNH4)
    if fluxdiff > 0 and prodiff > 0:
        zid = rxnid
    if fluxdiff < 0 and prodiff < 0:
        fid = rxnid
    return zid, fid, genename


def pro_flux(model, Nsource, flux, pro, SGD):
    zid_l = []
    fid_l = []
    zgene_l = []
    for r in model.reactions:
        gr = r.gene_reaction_rule
        a_s = gr.split(' and ')
        o_s = gr.split(' or ')
        rxnid = r.id
        if len(a_s) == 1 and len(o_s) == 1:
            # one to one
            zid, fid, pro_Nsource, pro_NH4, genename = onetoone(rxnid, gr, Nsource, flux, pro)
            if zid != 'nan':
                zid_l.append(zid)
                zgene_l.append(genename)
            if fid != 'nan':
                fid_l.append(fid)
        elif len(a_s) == 1 and len(o_s) != 1:
            # isoenzyme, or
            zid, fid, genename = isoenzyme(o_s, rxnid, Nsource, flux, pro, SGD)
            if zid != 'nan':
                zid_l.append(zid)
                for gii in genename:
                    zgene_l.append(gii)
            if fid != 'nan':
                fid_l.append(fid)
        elif len(a_s) != 1 and len(o_s) == 1:
            # complex
            zid, fid, pro_Nsource, pro_NH4, genename = complex(gr, rxnid, Nsource, flux, pro, SGD)
            if zid != 'nan':
                zid_l.append(zid)
                for gii in genename:
                    zgene_l.append(gii)
            if fid != 'nan':
                fid_l.append(fid)
        elif len(a_s) != 1 and len(o_s) != 1:
            # isoenzyme and complex
            zid, fid, genename = comiso(o_s, rxnid, Nsource, flux, pro, SGD)
            if zid != 'nan':
                zid_l.append(zid)
                for gii in genename:
                    zgene_l.append(gii)
            if fid != 'nan':
                fid_l.append(fid)

    return zid_l, fid_l, zgene_l


def get_flux():
    phe = pd.read_csv('../output/pfba_2111_phe.csv')
    ile = pd.read_csv('../output/pfba_2111_ile.csv')
    gln = pd.read_csv('../output/pfba_2111_gln.csv')
    nh4 = pd.read_csv('../output/pfba_2111_nh4.csv')
    flux = pd.DataFrame(index=phe.iloc[:, 0], columns=['Glu', 'NH4', 'Ile', 'Phe'])
    flux.loc[:, 'Glu'] = gln.iloc[:, 1].tolist()
    flux.loc[:, 'NH4'] = nh4.iloc[:, 1].tolist()
    flux.loc[:, 'Ile'] = ile.iloc[:, 1].tolist()
    flux.loc[:, 'Phe'] = phe.iloc[:, 1].tolist()
    return flux
# prepare figure data
flux = get_flux()
#flux.index = flux.loc[:, 'id']
pro = pd.read_csv('../data/proteomic.csv')
pro.index = pro.loc[:, 'Gene']
SGD = pd.read_table('../data/SGDgeneNames.tsv')
model = cobra.io.read_sbml_model('../data/yeast-GEM.xml')
SGD.index = SGD.loc[:, 'Systematic_name']
# glu_zid_l, glu_fid_l, glu_zgene_l = pro_flux(model, 'Glu', flux, pro, SGD)
phe_zid_l, phe_fid_l, phe_zgene_l = pro_flux(model, 'Phe', flux, pro, SGD)
ile_zid_l, ile_fid_l, ile_zgene_l = pro_flux(model, 'Ile', flux, pro, SGD)
phe_zgene_l = list(set(phe_zgene_l))
ile_zgene_l = list(set(ile_zgene_l))
phe_zgene_s = pd.Series(index=phe_zgene_l)
for p in phe_zgene_l:
    try:
        phe_zgene_s.loc[p] = SGD.loc[p, 'Standard_name']
    except:
        phe_zgene_s.loc[p] = p
ile_zgene_s = pd.Series(index=ile_zgene_l)
for il in ile_zgene_l:
    try:
        ile_zgene_s.loc[il] = SGD.loc[il, 'Standard_name']
    except:
        ile_zgene_s.loc[il] = il
phe_zgene_s.to_excel('../output/phe_zgene.xlsx', index=False)
ile_zgene_s.to_excel('../output/ile_zgene.xlsx', index=False)

same = [i for i in phe_zgene_l if i in ile_zgene_l]
pd.Series(same).to_excel('../output/Phe_Ile_same.xlsx', index=False)

# sub = pd.read_table('../data/uniqueSubsystems.tsv')
# sub.index = sub.loc[:, 'ID']
# phe = pd.DataFrame(index=phe_zid_l, columns=['subsystem'])
# for i in phe.index.values.tolist():
#     phe.loc[i, 'subsystem'] = sub.loc[i, 'subsystem_unique_@2021_12']
# phecount = phe.value_counts('subsystem')
# phefigure = pd.DataFrame(index=range(len(phecount.index)), columns=['subsystem', 'number'])
# phefigure.loc[:, 'subsystem'] = phecount.index
# phefigure.loc[:, 'number'] = phecount.loc[:].values.tolist()
# phe.to_excel('../output/phe_sub.xlsx')
#
# ile = pd.DataFrame(index=ile_zid_l, columns=['subsystem'])
# for i in ile.index.values.tolist():
#     ile.loc[i, 'subsystem'] = sub.loc[i, 'subsystem_unique_@2021_12']
# ilecount = ile.value_counts('subsystem')
# ilefigure = pd.DataFrame(index=range(len(ilecount.index)), columns=['subsystem', 'number'])
# ilefigure.loc[:, 'subsystem'] = ilecount.index
# ilefigure.loc[:, 'number'] = ilecount.loc[:].values.tolist()
# ile.to_excel('../output/ile_sub.xlsx')
#############################################
# figure
# font3 = {'family': 'Arial',
#          'weight': 'normal',
#          'size': 15,
#          }
# plt.figure(figsize=(20, 10))
# ax = plt.subplot(111, polar=True)
# plt.axis('off')
#
# upperLimit = 150
# lowerLimit = 5
# labelPadding = 4
#
# max = phefigure['number'].max()
#
# slope = (max - lowerLimit) / max
# heights = slope * phefigure.number + lowerLimit
#
# width = 2*np.pi / len(phefigure.index)
#
# indexes = list(range(1, len(phefigure.index)+1))
# angles = [element * width for element in indexes]
# angles
#
# # Draw bars
# bars = ax.bar(
#     x=angles,
#     height=heights,
#     width=width,
#     bottom=lowerLimit,
#     linewidth=2,
#     edgecolor="white",
#     color="#61a4b2",
# )
#
# # Add labels
# for bar, angle, height, label in zip(bars,angles, heights, phefigure["subsystem"]):
#
#     rotation = np.rad2deg(angle)
#
#     alignment = ""
#     if angle >= np.pi/2 and angle < 3*np.pi/2:
#         alignment = "right"
#         rotation = rotation + 180
#     else:
#         alignment = "left"
#
#     ax.text(
#         x=angle,
#         y=lowerLimit + bar.get_height() + labelPadding,
#         s=label,
#         ha=alignment,
#         va='center',
#         rotation=rotation,
#         rotation_mode="anchor",
#         fontdict=font3)
#
# plt.tight_layout()
# plt.savefig('../N_lim/output/phe.jpg',
#             dpi=600)
# plt.show()
#
# ###################################
# # ile
# font3 = {'family': 'Arial',
#          'weight': 'normal',
#          'size': 15,
#          }
# plt.figure(figsize=(20,10))
# ax = plt.subplot(111, polar=True)
# plt.axis('off')
#
# upperLimit = 150
# lowerLimit = 5
# labelPadding = 4
#
# max = ilefigure['number'].max()
#
# slope = (max - lowerLimit) / max
# heights = slope * ilefigure.number + lowerLimit
#
# width = 2*np.pi / len(ilefigure.index)
#
# indexes = list(range(1, len(ilefigure.index)+1))
# angles = [element * width for element in indexes]
# angles
#
# # Draw bars
# bars = ax.bar(
#     x=angles,
#     height=heights,
#     width=width,
#     bottom=lowerLimit,
#     linewidth=2,
#     edgecolor="white",
#     color="#61a4b2",
# )
#
# # Add labels
# for bar, angle, height, label in zip(bars,angles, heights, ilefigure["subsystem"]):
#
#     rotation = np.rad2deg(angle)
#
#     alignment = ""
#     if angle >= np.pi/2 and angle < 3*np.pi/2:
#         alignment = "right"
#         rotation = rotation + 180
#     else:
#         alignment = "left"
#
#     ax.text(
#         x=angle,
#         y=lowerLimit + bar.get_height() + labelPadding,
#         s=label,
#         ha=alignment,
#         va='center',
#         rotation=rotation,
#         rotation_mode="anchor",
#         fontdict=font3)
#
# plt.tight_layout()
# plt.savefig('../N_lim/output/ile.jpg',
#             dpi=600)
# plt.show()
