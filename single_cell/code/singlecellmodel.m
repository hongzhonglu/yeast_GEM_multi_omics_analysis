% get single cell model
% Data: Audrey P. Gasch 2017
% method: GIMME

model = importModel('../data/yeast-GEM.xml');
[num,txt] = xlsread('../data/GSE102475_GASCH_NaCl-scRNAseq_NormData.xlsx');
name = txt(1, (2:end));
changeCobraSolver ('glpk', 'all');
exp.gene = txt((2:end), 1);
model.lb(find(strcmp('r_0446', model.rxns()))) = -1000;
for i = 1 : length(name)
    exp.rawValue = num(:,i);
    exp.value = num(:,i);
    parsedGPR = get_GPR(model);
    [gene_id, gene_expr] = findUsedGenesLevels(model,exp);
    expressionRxns = selectGeneFromGPR(model, gene_id, gene_expr, parsedGPR, false);
    options.expressionRxns = expressionRxns;
    options.threshold = 20;
    options.solver = 'GIMME';
    cs_model = createTissueSpecificModel(model, options);
    modelFileName = convertStringsToChars(strcat('../output/scmodel/', name{i}, '.xml'));
    disp(output.flux{i})
    writeCbModel(cs_model, 'sbml', modelFileName);
    disp(i);
end

