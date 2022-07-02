%%  get pos_genes_in_react_expr
function [pos_genes_in_react_expr] = getPos(model)
fid = fopen('../data/yeast-GEM.txt');
GPR = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
RxnExp.rxnid = GPR{1};
RxnExp.gpr = GPR{3};
fclose(fid);
genes = model.genes;

for g = 1:length(genes)
    co = 1;
    for i = 1:length(RxnExp.gpr)
        tempgpr = split(RxnExp.gpr{i},[" and ", " or ", "(", ")"]);
        if any(strcmp(genes{g}, tempgpr), "all")
            temppos(co,1) = i;
            co = co + 1;
        end
    end

    pos_genes_in_react_expr{g} = temppos;
    clear temppos;
end
end

