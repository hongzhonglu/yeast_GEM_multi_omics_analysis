function parsedGPR = get_GPR(yeastmodel)
    parsedGPR = {};
    for i = 1:numel(yeastmodel.rxns)
        if ~isempty(yeastmodel.grRules{i})
            genes = split(yeastmodel.grRules{i}, ' or ');
            removeBracket = regexprep(genes,'[\(\)]', '');
            for j = 1:numel(removeBracket)
                gene = split(removeBracket{j}, ' and ');
                parsedGPR{i}{j}=gene;
            end
        else
            parsedGPR{i}={''};
        end
    end
    parsedGPR=parsedGPR';
end

