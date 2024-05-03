%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predict synthetic lethal
% Data: Michael Costanzo et al., Science(2016)

cd ../data/
% initCobraToolbox
disp('prepare data and model')
nxn = xlsread("data2016_in_model.xlsx");
generow = readmatrix("genelist_row.csv", 'OutputType', 'string');
genecol = readmatrix("genelist_col.csv", 'OutputType', 'string');
model = importModel("yeast-GEM.xml");
model = ravenCobraWrapper(model);
tp = 0; fp = 0; fn = 0; tn = 0;
cd ../code/
tic
doubleKO = {};
n = 1;
[~, grRateKO, ~] = doubleGeneDeletion(model,'FBA', generow, genecol);
for row = 1:length(generow)
    fprintf('%.4f \n', row/length(generow))
    for col = 1:length(genecol)
        growth = grRateKO(row, col);
            [tp, tn, fp, fn] = simu_exp(growth, nxn(row, col), tp, tn, fp, fn);
    end
end
fprintf('tp=%u\nfp=%u\ntn=%u\nfn=%u\n', tp, fp, tn, fn)
%writecell(doubleKO, '../output/2016genepair_fp_fn.xlsx')
disp('finish')
toc

