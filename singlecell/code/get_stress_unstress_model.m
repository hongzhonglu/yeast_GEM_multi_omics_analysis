modelPath = 'D:\model_research\yeast-GEM\model\yeast-GEM.xml';
model = readCbModel(modelPath);
[num,txt] = xlsread('D:\model_research\yeast_GEM_multi_omics_analysis\singlecell\data\stress_NaCl-scRNAseq_NormData.xlsx');
exp.gene = txt((2:end), 1);
exp.rawValue = num;
exp.value = mean(num, 2);
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, exp);
changeCobraSolver ('glpk', 'all');

options.expressionRxns = expressionRxns;
options.threshold = 50;
options.solver = 'GIMME';
GIMME_stress_all_model50 = createTissueSpecificModel(model, options);


