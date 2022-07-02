% open a yeast GEM
% run this code and get condition-specific flux
% the following data processing is in python
[glucosefluxes, glucosedead]=main(model,"glucose",'main');
disp('finish main');

[glucosefluxes, glucosedead]=main(model,"glucose", 'ITS');
disp('finish ITS');
