% open a yeast GEM
% run this code and get condition-specific flux
% the following data processing is in python
model = importModel('../data/yeast-GEM.xml');

[glucosefluxes, glucosedead]=main(model,"glucose",'main');
disp('finish main');

