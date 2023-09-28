import pandas as pd
import cobra
import logging


def change_medium(model, sub_type, sub):
    '''
        {'r_1714'}; % D-glucose exchange
        {'r_1654'}; % ammonium exchange
        {'r_2005'}; % phosphate exchange
        {'r_2060'}; % phosphate exchange
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
        logging.info('add {}'.format(rxn.reaction))
    except ValueError:
        met = model_test.metabolites.get_by_id(metid)
        print(metid)
        met.compartment = 'e'
        model_test.add_metabolites(met)
        rxn = model_test.add_boundary(met, type="exchange")
        logging.info('add {}'.format(rxn.reaction))
    return model_test, rxn

def simu(model):
    sol = model.optimize()
    if sol.objective_value > 1e-6:
        gng = 'G'
    else:
        gng = 'NG'
    return gng


logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')
biolog = pd.read_table('./Biolog_Substrate.tsv')
# yeast 9
# model = cobra.io.read_sbml_model('./yeast-GEM.xml')

# yeast8.3
model = cobra.io.read_sbml_model("./yeast-GEM-8.3.0/yeast-GEM-8.3.0/ModelFiles/xml/yeastGEM.xml")
exp_simu = pd.DataFrame(columns=['Name_in_Model', 'id', 'Substrate_type', 'Growth_Biolog', 'Growth_Model'])
for m in model.metabolites:
    for i in range(len(biolog.index)):
        #if m.name == biolog.loc[i, 'Name_in_Model']:
        if m.name.split(' ')[0] == biolog.loc[i, 'Name_in_Model']:
            exp_simu.loc[i, 'id'] = m.id
            exp_simu.loc[i, 'Substrate_type'] = biolog.loc[i, 'Substrate_type']
            exp_simu.loc[i, 'Growth_Biolog'] = biolog.loc[i, 'Growth_Biolog']
            exp_simu.loc[i, 'Name_in_Model'] = biolog.loc[i, 'Name_in_Model']

exp_simu.index = range(len(exp_simu))


if __name__ == '__main__':
    for i in range(len(exp_simu.index)):
        model_test, rxn = add_exchange(model, exp_simu.loc[i, 'id'])
        model_final = change_medium(model_test, exp_simu.loc[i, 'Substrate_type'], rxn)
        gng = simu(model_final)
        logging.info('{}: {}'.format(str(exp_simu.loc[i, 'Name_in_Model']), gng))
        exp_simu.loc[i, 'Growth_Model'] = gng

    c = 0
    for i in range(len(exp_simu.index)):
        if exp_simu.loc[i, 'Growth_Biolog'] == exp_simu.loc[i, 'Growth_Model']:
            c += 1

    exp_simu.to_excel('./{}_sub_use.xlsx'.format(model.id))
