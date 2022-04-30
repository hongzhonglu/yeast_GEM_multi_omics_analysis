# get flux data by simpling or pFBA
# from regulation_analysis_for_yeast9.get_fluxandproteoimcs.code.load_csmodel import load_csmodel
# [modelphe_Nlim_01, modelile_Nlim_01, modelNH4_Nlim_035, \
#            modelNH4_Nlim_030, modelNH4_Nlim_018, modelNH4_Nlim_013, \
#            modelCN115, modelCN50, modelCN30, modelNH4_Nlim_01, \
#            modelNH4_Nlim_005, modelgln_Nlim_01] = load_csmodel

import cobra
import pandas as pd


def getflux(model, analysis):
    # set ATPM as object
    model.objective = 'r_4046'
    model.reactions.get_by_id('r_4046').bounds = (0, 1000)
    solution = model.optimize()
    model.reactions.get_by_id('r_4046').bounds = (solution.objective_value, 1000)
    if analysis == 'sample':
        sam_solution = cobra.sampling.sample(model, 100)
        csflux = sam_solution.mean(axis=0)
    if analysis == 'pFBA':
        pfba_solution = cobra.flux_analysis.pfba(model)
        csflux = pfba_solution.fluxes
    return csflux


def getpro():
    CN_pro = pd.read_excel(r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\data\CNproteomics.xlsx')
    D_diff_pro = pd.read_csv(r'.\regulation_analysis_for_yeast9\get_fluxandproteoimcs\data\D_diff_proteomics.csv')
    return CN_pro, D_diff_pro
