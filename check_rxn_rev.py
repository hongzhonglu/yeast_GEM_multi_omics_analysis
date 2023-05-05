import cobra
import pandas as pd

model = cobra.io.read_sbml_model(r'D:\model_research\yeast-GEM\model\yeast-GEM.xml')
dgdata = pd.read_excel('D:\model_research\yeast_GEM_multi_omics_analysis\de_data_all.xlsx')
newmodel = model.copy()


for i in range(len(model.reactions)):
    if model.reactions[i].upper_bound == 1000.0 and model.reactions[i].lower_bound == -1000.0:
        if isinstance(dgdata.loc[i, 'detaG'], float) and dgdata.loc[i, 'detaG'] > 0:
            dgdata.loc[i, 'check'] = 'double check needed detaG > 0'
        elif isinstance(dgdata.loc[i, 'detaG'], float) and dgdata.loc[i, 'detaG'] < 0:
            dgdata.loc[i, 'check'] = 'double check needed detaG < -0'
        elif isinstance(dgdata.loc[i, 'detaG'], float):
            dgdata.loc[i, 'check'] = '0'
    elif model.reactions[i].upper_bound == 1000.0 and model.reactions[i].lower_bound == 0.0:
        if isinstance(dgdata.loc[i, 'detaG'], float) and dgdata.loc[i, 'detaG'] <= 0:
            dgdata.loc[i, 'check'] = 'right'
        elif isinstance(dgdata.loc[i, 'detaG'], float) and dgdata.loc[i, 'detaG'] > 0:
            dgdata.loc[i, 'check'] = 'model_F but detaG > 0'
    elif model.reactions[i].upper_bound == 0 and model.reactions[i].lower_bound == -1000:
        if isinstance(dgdata.loc[i, 'detaG'], float) and dgdata.loc[i, 'detaG'] < 0:
            dgdata.loc[i, 'check'] = 'model_R but detaG < 0'
        elif isinstance(dgdata.loc[i, 'detaG'], float) and dgdata.loc[i, 'detaG'] >= 0:
            dgdata.loc[i, 'check'] = 'right'

dgdata.to_excel('./dgcheck.xlsx')