# geneshort name is not contain in the RAVEN model
# so first run this python code to add that into the model.
def temp9():
    import cobra
    import pandas as pd
    #model = cobra.io.read_sbml_model(r"C:\Users\yuhuzhouye\Desktop\yeast-GEM.xml")
    model = cobra.io.read_sbml_model(r"D:\model_research\yeast-GEM\model\yeast-GEM.xml")
    SGD = pd.read_table(r'D:\model_research\yeast_GEM_multi_omics_analysis\synthetic lethal\data\SGDgeneNames.tsv')
    gene = model.genes
    SGD.index = SGD.loc[:, 'Systematic_name']
    for g in gene:
        s_name = SGD.loc[g.id, 'Standard_name']
        model.genes.get_by_id(g.id).name = s_name
    return model
