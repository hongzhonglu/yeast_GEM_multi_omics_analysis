%% The following codes are baesd on "www.pnas.org/cgi/doi/10.1073/pnas.2002959117"
function [fluxes, dead, model]=main(model,carbonsource, dataset)
[reaction_expression] = getRxnExp();
[pos_genes_in_react_expr] = getPos(model);
[idxs_genes_sorted_by_length] = getidxs(model);
[model] = ChangeCarbonSource(model,carbonsource);
% Select dataset ('main' or 'ITS' for the experimentally-independent test set)

% Load gene expression and growth data 
 if strcmp(dataset, 'main')
    load('MSB_data_expression.mat'); 
    load('MSB_data_growth_rates.mat');
    expr_data = msbdata;
    growth_rates = growthRatesMSB;
elseif strcmp(dataset, 'ITS')
    load('independent_test_set_expression.mat');
    load('independent_test_set_growth_rates.mat');
    expr_data = delmutantslimmaSameithetal2015;
    growth_rates = growthRatesSameithetal2015;
end

growthReactName = 'r_2111';
model = changeObjective(model, growthReactName); % Set the objetive to be the biomass
genes = model.genes;
genes_in_dataset = expr_data.geneName;
GeneExpressionArray = ones(numel(genes),1);
gamma = 1;
lowerPercentBounded = 1; % We initially explored different low bound values and settled on 1
bounds = -5;
model = setMediaConditions(model, bounds);
[model] = ChangeCarbonSource(model,carbonsource);
fluxes = model.rxns;
dead(1,1) = "lethal";
co = 1;

for t=3:(width(expr_data))
    expr_profile = table2array(expr_data(:,t));
    if contains(expr_data.Properties.VariableNames{t}, '_') % Double KO case
        deletedGene = strsplit(expr_data.Properties.VariableNames{t}, '_');
    else % Single KO case
        deletedGene = expr_data.Properties.VariableNames{t};
    end
    deletedGeneIndex = find(strcmp(expr_data.commonName,deletedGene),1);
    deletedGeneExpressionIndexInArray = inf;
    pos_genes_in_dataset = zeros(numel(genes),1);
    for i=1:numel(genes)
        position = find(strcmp(genes{i},genes_in_dataset));
        if ~isempty(position)
            pos_genes_in_dataset(i) = find(strcmp(genes{i},genes_in_dataset));
            if pos_genes_in_dataset(i) == deletedGeneIndex
                deletedGeneExpressionIndexInArray = i;
            end
            GeneExpressionArray(i) = expr_profile(pos_genes_in_dataset(i));
        end
    end

    NumberOfGenes = numel(genes);

    [v_out, f_out, LB, UB] = evaluate_objective(gamma,GeneExpressionArray,model,genes,reaction_expression,pos_genes_in_react_expr,idxs_genes_sorted_by_length,...
        lowerPercentBounded,deletedGeneExpressionIndexInArray);
    f_out

    if f_out > 0 
        fluxes = [fluxes, array2table(v_out,'VariableNames',{expr_data.Properties.VariableNames{t}})];
        modellb(:, co) = LB;
        modelub(:, co) = UB;
        co = co + 1;
    elseif f_out == 0
        dead(t-1,1) = string(expr_data.Properties.VariableNames{t});
    end
end
if strcmp(dataset, 'main')
    writetable(fluxes, strcat('../output/',carbonsource,'_fluxDataset.csv'), 'WriteRowNames', true);
    %writematrix(dead, strcat('../output/',carbonsource,'_deadDataset.csv'));
    %writematrix(modellb, strcat('../output/',carbonsource,'_LB.csv'));
    %writematrix(modelub, strcat('../output/',carbonsource,'_UB.csv'));
else
    writetable(fluxes, strcat('../output/Double',carbonsource,'_fluxDataset.csv'), 'WriteRowNames', true);
    %writematrix(dead, strcat('../output/Double',carbonsource,'_deadDataset.csv'));
    %writematrix(modellb, strcat('../output/Double',carbonsource,'_LB.csv'));
    %writematrix(modelub, strcat('../output/Double',carbonsource,'_UB.csv'));
end
 
end

