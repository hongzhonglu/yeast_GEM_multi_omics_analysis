import pandas as pd
import cobra


def change_medium(model, sub_type, sub):
    '''
        {'r_1714'}; % D-glucose exchange
        {'r_1654'}; % ammonium exchange
        {'r_2005'}; % phosphate exchange
        {'r_2060'}; % sulfur exchange
    '''
    model_test = model.copy()
    medium = model_test.medium
    medium[sub.id] = 10.0
    if sub_type == 'C':
        medium['r_1714'] = 0.0
    elif sub_type == 'N':
        medium['r_1654'] = 0.0
    elif sub_type == 'P':
        medium['r_2005'] = 0.0
    elif sub_type == 'S':
        medium['r_2060'] = 0.0
    model_test.medium = medium
    return model_test

def add_exchange(model, metid):
    model_test = model.copy()
    try:
        rxn = model_test.add_boundary(model_test.metabolites.get_by_id(metid))
    except ValueError:
        met = model_test.metabolites.get_by_id(metid)
        print(metid)
        met.compartment = 'e'
        model_test.add_metabolites(met)
        rxn = model_test.add_boundary(met, type="exchange")
    return model_test, rxn

def simu(model):
    try:
        sol = model.optimize()
    except UserWarning:
        # release ATPM constrains
        model.reactions.get_by_id('r_4046').lower_bound = 0.35
        sol = model.optimize()
    if sol.objective_value > 1e-6:
        gng = 'G'
    else:
        gng = 'NG'
    return gng


biolog = pd.read_table('./Biolog_Substrate.tsv')
# yeast 9
# model = cobra.io.read_sbml_model('./yeast-GEM.xml')

# yeast8.3
model = cobra.io.read_sbml_model("./yeast-GEM-8.3.0/yeast-GEM-8.3.0/ModelFiles/xml/yeastGEM.xml")
exp_simu = pd.DataFrame(columns=['Name_in_Model', 'id', 'Substrate_type', 'Growth_Biolog', 'Growth_Model'])
for m in model.metabolites:
    if m.compartment == 'c':
        for i in range(len(biolog.index)):
            # This "if" is for yeast8 because the metabolites name is like "name [compartment]"
            if m.name.split(' [')[0] == biolog.loc[i, 'Name_in_Model']:
            # This "if" is for yeast9
            # if m.name == biolog.loc[i, 'Name_in_Model']:
                exp_simu.loc[i, 'id'] = m.id
                exp_simu.loc[i, 'Substrate_type'] = biolog.loc[i, 'Substrate_type']
                exp_simu.loc[i, 'Growth_Biolog'] = biolog.loc[i, 'Growth_Biolog']
                exp_simu.loc[i, 'Name_in_Model'] = biolog.loc[i, 'Name_in_Model']

exp_simu.index = range(len(exp_simu))


if __name__ == '__main__':
    allrxnName = [r.name for r in model.reactions]
    allrxn = [r for r in model.reactions]
    for i in range(len(exp_simu.index)):
        if exp_simu.loc[i, 'Substrate_type'] != 'S':
            model_test, rxn = add_exchange(model, exp_simu.loc[i, 'id'])
            model_final = change_medium(model_test, exp_simu.loc[i, 'Substrate_type'], rxn)
            gng = simu(model_final)
            exp_simu.loc[i, 'Growth_Model'] = gng
        else:
            try:
                name = exp_simu.loc[i, 'Name_in_Model'] + ' exchange'
                rxn = allrxn[allrxnName.index(name)]
                model_final = change_medium(model, exp_simu.loc[i, 'Substrate_type'], rxn)
                gng = simu(model_final)
                exp_simu.loc[i, 'Growth_Model'] = gng
            except:
                pass
    c = 0
    for i in range(len(exp_simu.index)):
        if exp_simu.loc[i, 'Growth_Biolog'] == exp_simu.loc[i, 'Growth_Model']:
            c += 1

    exp_simu.to_excel('./{}_sub_use_SolveCompartmentIssue.xlsx'.format(model.id))
