% get single cell model
% Data: Audrey P. Gasch 2017
% method: GIMME

model = readCbModel('../data/yeast-GEM.xml');
[num,txt] = xlsread('../data/GSE102475_GASCH_NaCl-scRNAseq_NormData.xlsx');
name = txt(1, (2:end));
changeCobraSolver ('glpk', 'all');
exp.gene = txt((2:end), 1);
for i = 1 : length(name)
    exp.rawValue = num(:,i);
    exp.value = num(:,i);
    [expressionRxns, parsedGPR] = mapExpressionToReactions(model, exp);
    options.expressionRxns = expressionRxns;
    options.threshold = 24;
    options.solver = 'GIMME';
    cs_model = createTissueSpecificModel(model, options);
    modelFileName = convertStringsToChars(strcat('../output/scmodel/', name{i}, '.xml'));
    writeCbModel(cs_model, 'sbml', modelFileName);
    disp(i);
end
 
