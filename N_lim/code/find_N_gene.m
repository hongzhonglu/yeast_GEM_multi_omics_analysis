model = importModel("C:\Users\yuhuzhouye\Desktop\yeast9_w\yesat9_temp\yeast-GEM.xml");
target_sub = {
    'Nitrogen metabolism';
    'Purine metabolism';
    'Pyrimidine metabolism';
    'Alanine, aspartate and glutamate metabolism';
    'Glycine, serine and threonine metabolism';
    'Cysteine and methionine metabolism';
    'Valine, leucine and isoleucine degradation';
    'Valine, leucine and isoleucine biosynthesis';
    'Lysine biosynthesis';
    'Lysine degradation';
    'Arginine biosynthesis';
    'Arginine and proline metabolism';
    'Histidine metabolism';
    'Tyrosine metabolism';
    'Phenylalanine metabolism';
    'Tryptophan metabolism';
    'Phenylalanine, tyrosine and tryptophan biosynthesis'
    };
rxn = {};
rxn{1, 1} = 'rxn_name';
rxn{1, 2} = 'gpr';
rxn{1, 3} = 'sub';
c = 1;
for i = 1:length(model.subSystems)
    if ~isempty(model.grRules{i, 1}) && ~isempty(find(strcmp(cell2str(model.subSystems{i}), target_sub), 1))
        rxn{c+1, 1} = model.rxnNames{i, 1};
        rxn{c+1, 2} = model.grRules{i, 1};
        rxn{c+1, 3} = cell2str(model.subSystems{i, 1});
        c = c + 1;
    end
end
writecell(rxn, '../output/N_rxns.csv')
