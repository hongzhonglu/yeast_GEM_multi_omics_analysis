%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predict synthetic lethal
% Data: Michael Costanzo et al., Science(2017)

cd ../data/
% initCobraToolbox
disp('prepare data and model')
nxn = xlsread("data2016_in_model.xlsx");
generow = readmatrix("genelist_row.csv", 'OutputType', 'string');
genecol = readmatrix("genelist_col.csv", 'OutputType', 'string');
model = importModel("yeast-GEM.xml");

cd ../code/
[~, grRateKO, ~] = doubleGeneDeletion(model,'FBA', generow, genecol);
for row = 1:length(generow)
    for col = 1:length(genecol)
         [tp, tn, fp, fn] = simu_exp(grRateKO(row, col), nxn(row, col), tp, tn, fp, fn);
    end
end
fprintf('tp=%u\nfp=%u\ntn=%u\nfn=%u\n', tp, fp, tn, fn)
disp('finish')


