# get flux data by simpling or pFBA
# read proteomics data
# [modelphe_Nlim_01, modelile_Nlim_01, modelNH4_Nlim_035, \
#            modelNH4_Nlim_030, modelNH4_Nlim_018, modelNH4_Nlim_013, \
#            modelCN115, modelCN50, modelCN30, modelNH4_Nlim_01, \
#            modelNH4_Nlim_005, modelgln_Nlim_01] = load_csmodel

import cobra
import pandas as pd
from regulation_analysis_for_yeast9.get_fluxandproteoimcs.code.load_csmodel import load_csmodel


def getflux(model, analysis):
    # model   :   metabolic model
    # analysis:   pFBA or sample
    # set ATPM as object
    model.objective = 'r_4046'
    model.reactions.get_by_id('r_4046').bounds = (0, 1000)
    solution = model.optimize()
    model.reactions.get_by_id('r_4046').bounds = (solution.objective_value, 1000)
    if analysis == 'sample':
        sam_solution = cobra.sampling.sample(model, 100)
        csflux = sam_solution.mean(axis=0).to_dict()
        gene_pro = []
        for r in csflux.keys():
            gene_pro.append(model.reactions.get_by_id(r).gene_name_reaction_rule)
        flux_gene_value = pd.DataFrame({
            'reactionID': csflux.keys(),
            'flux': csflux.values(),
            'gene': gene_pro,
        })
    if analysis == 'pFBA':
        pfba_solution = cobra.flux_analysis.pfba(model)
        csflux = pfba_solution.fluxes.to_dict()
        gene_pro = []
        for r in csflux.keys():
            gene_pro.append(model.reactions.get_by_id(r).gene_name_reaction_rule)
        flux_gene_value = pd.DataFrame({
            'reactionID':csflux.keys(),
            'flux': csflux.values(),
            'gene': gene_pro},
            )
    flux_gene_value.index = flux_gene_value.loc[:, 'reactionID']
    return flux_gene_value


def getpro(condition):
    # get proteomics data under specific condition
    # condition = CN30, CN50, CN115, N005, N010, N013, N018, N030, N035, Gln, Phe, Ile
    CN_pro = pd.read_excel(r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\data\CNproteomics.xlsx')
    D_diff_pro = pd.read_csv(r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\data\D_diff_proteomics.csv')
    if condition == 'CN30':
        pro_means = []
        for p in CN_pro.index:
            temp = (CN_pro.loc[p, 'C/N=30 rep1 prot (fmol mgDW-1)'] +
                    CN_pro.loc[p, 'C/N=30 rep2 prot (fmol mgDW-1)']) / 2
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': CN_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'CN50':
        pro_means = []
        for p in CN_pro.index:
            temp = (CN_pro.loc[p, 'C/N=50 rep1 prot (fmol mgDW-1)'] +
                    CN_pro.loc[p, 'C/N=50 rep2 prot (fmol mgDW-1)']) / 2
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': CN_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'CN115':
        pro_means = []
        for p in CN_pro.index:
            temp = (CN_pro.loc[p, 'C/N=115 rep1 prot (fmol mgDW-1)'] +
                    CN_pro.loc[p, 'C/N=115 rep2 prot (fmol mgDW-1)']) / 2
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': CN_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'N005':
        pro_means = []
        for p in D_diff_pro.index:
            temp = (D_diff_pro.loc[p, 'prot.1'] +
                    D_diff_pro.loc[p, 'prot.2'] +
                    D_diff_pro.loc[p, 'prot.3']
                    ) / 3
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': D_diff_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'N010':
        pro_means = []
        for p in D_diff_pro.index:
            temp = (D_diff_pro.loc[p, 'prot.7'] +
                    D_diff_pro.loc[p, 'prot.8'] +
                    D_diff_pro.loc[p, 'prot.9']
                    ) / 3
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': D_diff_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'N013':
        pro_means = []
        for p in D_diff_pro.index:
            temp = (D_diff_pro.loc[p, 'prot.10'] +
                    D_diff_pro.loc[p, 'prot.11'] +
                    D_diff_pro.loc[p, 'prot.12']
                    ) / 3
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': D_diff_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'N018':
        pro_means = []
        for p in D_diff_pro.index:
            temp = (D_diff_pro.loc[p, 'prot.13'] +
                    D_diff_pro.loc[p, 'prot.14'] +
                    D_diff_pro.loc[p, 'prot.15']
                    ) / 3
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': D_diff_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'N030':
        pro_means = []
        for p in D_diff_pro.index:
            temp = (D_diff_pro.loc[p, 'prot.16'] +
                    D_diff_pro.loc[p, 'prot.17'] +
                    D_diff_pro.loc[p, 'prot.18']
                    ) / 3
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': D_diff_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'N035':
        pro_means = []
        for p in D_diff_pro.index:
            temp = (D_diff_pro.loc[p, 'prot.19'] +
                    D_diff_pro.loc[p, 'prot.20'] +
                    D_diff_pro.loc[p, 'prot.21']
                    ) / 3
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': D_diff_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'Gln':
        pro_means = []
        for p in D_diff_pro.index:
            temp = (D_diff_pro.loc[p, 'prot.28'] +
                    D_diff_pro.loc[p, 'prot.29'] +
                    D_diff_pro.loc[p, 'prot.30']
                    ) / 3
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': D_diff_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'Phe':
        pro_means = []
        for p in D_diff_pro.index:
            temp = (D_diff_pro.loc[p, 'prot.34'] +
                    D_diff_pro.loc[p, 'prot.35'] +
                    D_diff_pro.loc[p, 'prot.36']
                    ) / 3
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': D_diff_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    if condition == 'Ile':
        pro_means = []
        for p in D_diff_pro.index:
            temp = (D_diff_pro.loc[p, 'prot.40'] +
                    D_diff_pro.loc[p, 'prot.41'] +
                    D_diff_pro.loc[p, 'prot.42']
                    ) / 3
            pro_means.append(temp)
        pro_gene_values = pd.DataFrame({
            'gene': D_diff_pro.loc[:, 'Gene'],
            'prot': pro_means,
        })
    return pro_gene_values
