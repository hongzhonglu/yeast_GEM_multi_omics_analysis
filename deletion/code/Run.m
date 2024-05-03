% open a yeast GEM
% run this code and get condition-specific flux
model = importModel('../data/yeast-GEM.xml');

[glucosefluxes, glucosedead]=main(model,"glucose",'main');
[glucosefluxes, glucosedead]=main(model,"glucose",'ITS');
disp('finish main');

