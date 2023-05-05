model = loadYeastModel;
m03 = readcell('..\output\m05.csv');
for i = 2:height(m05)
    m03(i, 3) = model.subSystems(find(strcmp(m03(i, 2), model.rxns)));
end

tca = {};
for i = length(model.rxns)
    tca = model.rxns(find(strcmp('Citrate cycle (TCA cycle)', model.subSystems)));
end
a = cell2table(model.subSystems);

writetable(a,"D:\model_research\yeast_GEM_multi_omics_analysis\N_lim\data\yeast9_with_sub.xlsx")
