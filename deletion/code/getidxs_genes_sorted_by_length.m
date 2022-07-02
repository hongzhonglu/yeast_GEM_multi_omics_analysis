model = loadYeastModel;
genes = model.genes;
n5 = 1;
n7 = 1;
n9 = 1;
for i = 1:length(genes)
    if length(genes{i}) == 5
        g5{n5} = genes{i};
        n5 = n5 + 1;
    elseif length(genes{i}) == 7
        g7{n7} = genes{i};
        n7 = n7 + 1;
    elseif length(genes{i}) == 9
        g9{n9} = genes{i};
        n9 = n9 + 1;
    end
end

disp(length(g9)+length(g7)+length(g5)==length(genes))

sort_g = [g9, g7, g5];

for j=1:length(sort_g)
    idxs_genes_sorted_by_length(1,j) = find(strcmp(sort_g{j}, genes));
end