import cobra


def simu(model, source, method):
    if method == 'sample':
        print(model.id)
        model.objective = 'r_2111'
        model.reactions.get_by_id('r_2111').bounds = (0, 1000)
        solution = cobra.sampling.sample(model, 1000)
        solution.to_csv('../output/fig4_Sampling_{}.csv'.format(source))
        meanflux = solution.mean(axis=0)
        meanflux.to_csv('../output/fig4_Sampling_mean_{}.csv'.format(source))
    if method =='pFBA':
        print(model.id)
        model.objective = 'r_2111'
        solution = cobra.flux_analysis.pfba(model)
        solution.fluxes.to_csv('../output/pfba_2111_{}.csv'.format(source))
        meanflux = ''
    return solution, meanflux


modelgln_Nlim_01 = cobra.io.read_sbml_model(
    r'../output/Nlim_model/modelgln_Nlim_01.xml')
modelile_Nlim_01 = cobra.io.read_sbml_model(
    r'../output/Nlim_model/modelile_Nlim_01.xml')
modelNH4_Nlim_01 = cobra.io.read_sbml_model(
    r'../output/Nlim_model/modelNH4_Nlim_01.xml')
modelphe_Nlim_01 = cobra.io.read_sbml_model(
    r'../output/Nlim_model/modelphe_Nlim_01.xml')


sam_solution, meanflux = simu(modelgln_Nlim_01, 'gln', 'pFBA')
sam_solution1, meanflux1 = simu(modelile_Nlim_01, 'ile', 'pFBA')
sam_solution2, meanflux2 = simu(modelNH4_Nlim_01, 'nh4', 'pFBA')
sam_solution3, meanflux3 = simu(modelphe_Nlim_01, 'phe', 'pFBA')
