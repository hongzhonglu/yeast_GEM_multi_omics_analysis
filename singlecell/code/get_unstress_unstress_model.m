modelPath = 'D:\model_research\yeast-GEM\model\yeast-GEM.xml';
model = readCbModel(modelPath);
[num,txt] = xlsread('D:\model_research\yeast_GEM_multi_omics_analysis\singlecell\data\GSE102475_GASCH_NaCl-scRNAseq_NormData.xlsx');
name = txt(1, (2:end));
changeCobraSolver ('glpk', 'all');
exp.gene = txt((2:end), 1);

for i = 1 to length(name)
    exp.rawValue = num(:,i);
    exp.value = num(:,i);
    [expressionRxns, parsedGPR] = mapExpressionToReactions(model, exp);
end

options.expressionRxns = expressionRxns;
options.threshold = 24;
options.solver = 'GIMME';
GIMME_unstress_all_model3 = createTissueSpecificModel(model, options);


