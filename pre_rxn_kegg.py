# get rxn with full kegg annotation for dGpredictor

import cobra
import pandas as pd
import re



def get_eq(model, equation):
    ss = re.compile(r's_\d{4}')
    sss = ss.findall(equation)
    for s in sss:
        try:
            equation = equation.replace(s, model.metabolites.get_by_id(s).annotation['kegg.compound'])
            equation = equation.replace('-->', '<=>')
        except:
            equation = 'nan'
    return equation


model = cobra.io.read_sbml_model(r"D:\model_research\yeast-GEM\model\yeast-GEM.xml")
tmprxnEquation = [e.reaction for e in model.reactions]
rxnEquation = [get_eq(model, eq) for eq in tmprxnEquation]

dgpredictor_rxn = pd.DataFrame()
dgpredictor_rxn.loc[:, 'eq'] = rxnEquation
dgpredictor_rxn.to_excel('./dgpredictor.xlsx', index=False)


