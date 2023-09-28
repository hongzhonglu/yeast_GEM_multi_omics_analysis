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
in_model_row = {};
[~, grRateKO, ~] = doubleGeneDeletion(model,'FBA', generow, genecol);
dead_simu = find(grRateKO <= 0.00001);
dead_exp = find(-1000 < nxn <= -0.35);
if length(dead_simu) == length(dead_exp)
    disp('all true')
end



