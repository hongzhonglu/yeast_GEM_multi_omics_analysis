import cobra
import pandas as pd

def double_deletion(model, gene1, gene2, opt):
    modeltemp = model.copy()
    for r in model.reactions:
        genesplit = r.gene_reaction_rule.split(' ')
        if gene1 in genesplit \
                and gene2 not in genesplit \
                and 'or' not in genesplit:
            modeltemp.reactions.get_by_id(r.id).bounds = (0.0, 0.0)

        if gene2 in genesplit \
                and gene1 not in genesplit \
                and 'or' not in genesplit:
            modeltemp.reactions.get_by_id(r.id).bounds = (0.0, 0.0)

        if gene2 in genesplit \
                and gene1 in genesplit \
                and 'or' not in genesplit:
            modeltemp.reactions.get_by_id(r.id).bounds = (0.0, 0.0)

        if gene2 in genesplit \
                and gene1 in genesplit \
                and 'or' in genesplit \
                and len(genesplit) > 5:
            modeltemp.reactions.get_by_id(r.id).bounds = (0.0, 0.0)

    solu = modeltemp.optimize()
    if solu.objective_value < opt:
        live = 'negative'
    elif solu.objective_value > opt:
        live = 'positive'
    else:
        live = 'nan'
    return live

model = cobra.io.read_sbml_model(r'synthetic lethal/data/yeast-GEM.xml')
exp = pd.read_table(r'synthetic lethal/data/SGA_NxN_clustered.tsv')
exp.index = exp.loc[:, 'A']
lie = exp.index.values.tolist()
hang = exp.columns.values.tolist()
modelgene = [g.name for g in model.genes]

pair1 = [p1 for p1 in lie if p1 in modelgene]
pair2 = [p2 for p2 in hang if p2 in modelgene]

simu = []
ex = []
opt = model.optimize()
for pp1 in pair1[5:10]:
    for pp2 in pair2:
        live = double_deletion(model, pp1, pp2, opt.objective_value)
        simu.append(live)
        if exp.loc[pp1, pp2] < 0:
            ex.append('negative')
        elif exp.loc[pp1, pp2] > 0:
            ex.append('positive')
        else:
            ex.append('nan')
