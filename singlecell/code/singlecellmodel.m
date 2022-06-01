modelPath = 'D:\model_research\yeast-GEM\model\yeast-GEM.xml';
model = readCbModel(modelPath);
[num,txt] = xlsread('D:\model_research\yeast_GEM_multi_omics_analysis\singlecell\data\GSE102475_GASCH_NaCl-scRNAseq_NormData.xlsx');
name = txt(1, (2:end));
changeCobraSolver ('glpk', 'all');
exp.gene = txt((2:end), 1);
for i = 161 : length(name)
    exp.rawValue = num(:,i);
    exp.value = num(:,i);
    [expressionRxns, parsedGPR] = mapExpressionToReactions(model, exp);
    options.expressionRxns = expressionRxns;
    options.threshold = 24;
    options.solver = 'GIMME';
    cs_model = createTissueSpecificModel(model, options);
    modelFileName = convertStringsToChars(strcat('D:\model_research\yeast_GEM_multi_omics_analysis\singlecell\output\', name{i}, '.xml'));
    writeCbModel(cs_model, 'sbml', modelFileName);
end
 
