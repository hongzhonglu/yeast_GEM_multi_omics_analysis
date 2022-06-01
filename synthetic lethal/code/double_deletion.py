
def double_deletion(model, gene1, gene2):
    for r in model.reactions:
        genesplit = r.gene_reaction_rule.split(' ')
        if gene1 in genesplit \
                and gene2 not in genesplit \
                and 'or' not in genesplit:
            model.reactions.get_by_id(r.id).bounds = (0.0, 0.0)

        if gene2 in genesplit \
                and gene1 not in genesplit \
                and 'or' not in genesplit:
            model.reactions.get_by_id(r.id).bounds = (0.0, 0.0)

        if gene2 in genesplit \
                and gene1 in genesplit \
                and 'or' not in genesplit:
            model.reactions.get_by_id(r.id).bounds = (0.0, 0.0)

        if gene2 in genesplit \
                and gene1 in genesplit \
                and 'or' in genesplit \
                and len(genesplit) > 5:
            model.reactions.get_by_id(r.id).bounds = (0.0, 0.0)

    solu = model.optimize()
    if solu.objective_value < 0.000001:
        live = 'SL'
    else:
        live = 'non_SL'
    return live
