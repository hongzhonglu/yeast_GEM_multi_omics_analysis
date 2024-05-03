%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run double gene deletion
% gene is in systematic name
function [growth, g1, g2] = doubleGD(model, gene1, gene2)
    ind1 = find(strcmp(gene1, model.genes));
    ind2 = find(strcmp(gene2, model.genes));
    if ~isempty(ind1) && ~isempty(ind2)
        %fprintf('deleting %s and %s\n', gene1, gene2)
        [~, grRateKO, ~] = doubleGeneDeletion(model,'FBA', model.genes(ind1), model.genes(ind2));
        growth = grRateKO;
        g1 = gene1;
        g2 = gene2;
        %disp(growth)
    else
        %fprintf('%s or %s not in model\n', gene1, gene2)
        growth = 'nan';
        g1 = "nan";
        g2 = "nan";
    end
end









