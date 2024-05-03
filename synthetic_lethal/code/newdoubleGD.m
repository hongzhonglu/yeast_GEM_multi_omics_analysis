%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code handles double gene deletion regradless of "AND/OR"
% relationship, following reviewer's suggestions.
function [growth, g1, g2] = newdoubleGD(model, gene1, gene2)
ind1 = find(strcmp(gene1, model.genes));
ind2 = find(strcmp(gene2, model.genes));
if ~isempty(ind1) && ~isempty(ind2)
    for i = 1: length(model.rxns)
        if strfind(model.grRules{i}, gene1) ~= 0
            model.lb(i) = 0;
            model.ub(i) = 0;
            %fprintf('Delet rxn %s: %s \n', model.rxns{i}, gene1);
        end
        if strfind(model.grRules{i}, gene2) ~= 0
            model.lb(i) = 0;
            model.ub(i) = 0;
            %fprintf('Delet rxn %s: %s \n', model.rxns{i}, gene2);
        end
    end
    sol = optimizeCbModel(model);
    growth = sol.f;
    g1 = gene1;
    g2 = gene2;
else
    fprintf('%s or %s not in model\n', gene1, gene2)
    growth = 'nan';
    g1 = "nan";
    g2 = "nan";
end
end


