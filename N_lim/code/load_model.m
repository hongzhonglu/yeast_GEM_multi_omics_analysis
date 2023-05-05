% get N lim condition-specific model based on rosemarys' dataset
function [cs_model] = load_model
cd ..
model = importModel('data/yeast-GEM.xml');
exp1e = readmatrix('data/1e mean.xlsx');
si4 = readmatrix('data/Supplementary Data 4.xlsx');
modelNH4_Nlim_01 = model; modelgln_Nlim_01 = model; modelphe_Nlim_01 = model; modelile_Nlim_01 = model;
modelNH4_Nlim_005 = model; modelNH4_Nlim_013 = model; modelNH4_Nlim_018 = model; modelNH4_Nlim_030 =model; modelNH4_Nlim_035 = model;
modelCN30 = model; modelCN50 = model; modelCN115 = model;
% scaleBioMass is a function in yeast-GEM
% cd C:\Users\yuhuzhouye\Desktop\research-2022\com\yeast-GEM\code\otherChanges

%% NH4_Nlim_01
% constrain growth, exchange reaction, biomass
modelNH4_Nlim_01 = scaleBioMass(modelNH4_Nlim_01, 'protein', 0.22557);
modelNH4_Nlim_01 = scaleBioMass(modelNH4_Nlim_01, 'RNA', 0.0253);
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_2111',	0.1, 'b');
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1714',	-1.503, 'b');%-1.503
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1654',	-1000, 'l'); 
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1672',	exp1e(7,11),'b');%co2
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1992',	-exp1e(7,12),'b');%o2
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1761',	exp1e(7,13),'b');%ethanol
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1634',	exp1e(7,14),'b');%acetate 
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_2033',	exp1e(7,15),'b');%pyruvate 
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_2056',	exp1e(7,16),'b');%succinate 
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1808',	exp1e(7,17),'b');%glycerol 
modelNH4_Nlim_01 = changeObjective(modelNH4_Nlim_01,	'r_1654');
solutionNH4_Nlim_01 = optimizeCbModel(modelNH4_Nlim_01,'max');
% constrain the minimized N usage
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1654',	solutionNH4_Nlim_01.f, 'b');%0.39
cs_model.modelNH4_Nlim_01 = modelNH4_Nlim_01;

%% gln_Nlim_01
modelgln_Nlim_01 = scaleBioMass(modelgln_Nlim_01,'protein',0.30677);
modelgln_Nlim_01 = scaleBioMass(modelgln_Nlim_01,'RNA',0.02215);
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_2111',	0.1,	'b');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1714',	-1.313,	'b');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1654',	0,  'b');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1891',	-1000, 'l');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1672',	exp1e(9,11),'b');%co2
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1992',	-exp1e(9,12),'b');%o2
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1761',	exp1e(9,13),'b');%ethanol
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1634',	exp1e(9,14),'b');%acetate 
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_2033',	exp1e(9,15),'b');%pyruvate 
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_2056',	exp1e(9,16),'b');%succinate 
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1808',	exp1e(9,17),'b');%glycerol 
modelgln_Nlim_01 = changeObjective(modelgln_Nlim_01,	'r_1891');
solutiongln_Nlim_01 = optimizeCbModel(modelgln_Nlim_01,'max');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1891',	solutiongln_Nlim_01.f, 'b');%-0.56
cs_model.modelgln_Nlim_01 = modelgln_Nlim_01;

%% phe_Nlim_01
modelphe_Nlim_01 = scaleBioMass(modelphe_Nlim_01,'protein',0.3985);
modelphe_Nlim_01 = scaleBioMass(modelphe_Nlim_01,'RNA',0.02423);
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_2111',	0.1,	'b');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1714',	-1.673,	'b');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1654',	0,'b');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1903',	-1000,'l');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1672',	exp1e(11,11),'b');%co2
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1992',	-exp1e(11,12),'b');%o2
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1761',	exp1e(11,13),'b');%ethanol
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1634',	exp1e(11,14),'b');%acetate 
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_2033',	exp1e(11,15),'b');%pyruvate 
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_2056',	exp1e(11,16),'b');%succinate 
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1808',	exp1e(11,17),'b');%glycerol 
modelphe_Nlim_01 = changeObjective(modelphe_Nlim_01,	'r_1903');
solutionphe_Nlim_01 = optimizeCbModel(modelphe_Nlim_01,'max');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1903',	solutionphe_Nlim_01.f, 'b');%-0.59
cs_model.modelphe_Nlim_01 = modelphe_Nlim_01;

%% ile_Nlim_01
modelile_Nlim_01 = scaleBioMass(modelile_Nlim_01,'protein',0.43732);
modelile_Nlim_01 = scaleBioMass(modelile_Nlim_01,'RNA',0.02381);
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_2111',	0.1,	'b');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1714',	-1.594,	'b');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1654',	0,'b');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1897',	-1000,'l');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1672',	exp1e(13,11),'b');%co2
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1992',	-exp1e(13,12),'b');%o2
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1761',	exp1e(13,13),'b');%ethanol
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1634',	exp1e(13,14),'b');%acetate 
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_2033',	exp1e(13,15),'b');%pyruvate 
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_2056',	exp1e(13,16),'b');%succinate
% the tempYeast9 is infeasible when constrain the following exchange
% reaction
% modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1808',	exp1e(13,17),'b');%glycerol 
modelile_Nlim_01 = changeObjective(modelile_Nlim_01,	'r_1897');
solutionile_Nlim_01 = optimizeCbModel(modelile_Nlim_01,'max');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1897',	solutionile_Nlim_01.f, 'b');%-2.89
cs_model.modelile_Nlim_01 = modelile_Nlim_01;

%% NH4_Nlim_005
modelNH4_Nlim_005 = scaleBioMass(modelNH4_Nlim_005,'protein',0.26148);
modelNH4_Nlim_005 = scaleBioMass(modelNH4_Nlim_005,'RNA',0.02417);
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_2111',	0.05,	'b');
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1714',	-0.765,	'b');
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1654',	-1000,'l');
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1672',	exp1e(1,11),'b');%co2
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1992',	-exp1e(1,12),'b');%o2
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1761',	exp1e(1,13),'b');%ethanol
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1634',	exp1e(1,14),'b');%acetate 
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_2033',	exp1e(1,15),'b');%pyruvate 
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_2056',	exp1e(1,16),'b');%succinate 
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1808',	exp1e(1,17),'b');%glycerol 
modelNH4_Nlim_005 = changeObjective(modelNH4_Nlim_005,	'r_1654');
solutionNH4_Nlim_005 = optimizeCbModel(modelNH4_Nlim_005,'max');
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1654',	solutionNH4_Nlim_005.f, 'b');%-0.21
cs_model.modelNH4_Nlim_005 = modelNH4_Nlim_005;

%% NH4_Nlim_013
modelNH4_Nlim_013 = scaleBioMass(modelNH4_Nlim_013,'protein',0.27569);
modelNH4_Nlim_013 = scaleBioMass(modelNH4_Nlim_013,'RNA',0.02583);
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_2111',	0.13,	'b');
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1714',	-2.071,	'b');
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1654',	-1000,'l');
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1672',	exp1e(3,11),'b');%co2
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1992',	-exp1e(3,12),'b');%o2
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1761',	exp1e(3,13),'b');%ethanol
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1634',	exp1e(3,14),'b');%acetate 
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_2033',	exp1e(3,15),'b');%pyruvate 
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_2056',	exp1e(3,16),'b');%succinate 
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1808',	exp1e(3,17),'b');%glycerol 
modelNH4_Nlim_013 = changeObjective(modelNH4_Nlim_013,	'r_1654');
solutionNH4_Nlim_013 = optimizeCbModel(modelNH4_Nlim_013,'max');
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1654',	solutionNH4_Nlim_013.f, 'b');%-0.59
cs_model.modelNH4_Nlim_013 = modelNH4_Nlim_013;

%% NH4_Nlim_018
modelNH4_Nlim_018 = scaleBioMass(modelNH4_Nlim_018,'protein',0.28647);
modelNH4_Nlim_018 = scaleBioMass(modelNH4_Nlim_018,'RNA',0.03854);
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_2111',	0.18,	'b');
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1714',	-3.099,	'b');
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1654',	-1000,'l');
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1672',	exp1e(4,11),'b');%co2
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1992',	-exp1e(4,12),'b');%o2
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1761',	exp1e(4,13),'b');%ethanol
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1634',	exp1e(4,14),'b');%acetate 
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_2033',	exp1e(4,15),'b');%pyruvate 
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_2056',	exp1e(4,16),'b');%succinate 
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1808',	exp1e(4,17),'b');%glycerol 
modelNH4_Nlim_018 = changeObjective(modelNH4_Nlim_018,	'r_1654');
solutionNH4_Nlim_018 = optimizeCbModel(modelNH4_Nlim_018,'max');
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1654',	solutionNH4_Nlim_018.f, 'b');%-0.86
cs_model.modelNH4_Nlim_018 = modelNH4_Nlim_018;

%% NH4_Nlim_030
modelNH4_Nlim_030 = scaleBioMass(modelNH4_Nlim_030,'protein',0.41448);
modelNH4_Nlim_030 = scaleBioMass(modelNH4_Nlim_030,'RNA',0.04822);
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_2111',	0.30,	'b');
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1714',	-8.40,	'b');
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1654',	-1000,'l');
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1672',	exp1e(5,11),'b');%co2
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1992',	-exp1e(5,12),'b');%o2
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1761',	exp1e(5,13),'b');%ethanol
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1634',	exp1e(5,14),'b');%acetate 
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_2033',	exp1e(5,15),'b');%pyruvate 
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_2056',	exp1e(5,16),'b');%succinate 
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1808',	exp1e(5,17),'b');%glycerol 
modelNH4_Nlim_030 = changeObjective(modelNH4_Nlim_030,	'r_1654');
solutionNH4_Nlim_030 = optimizeCbModel(modelNH4_Nlim_030,'max');
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1654',	solutionNH4_Nlim_030.f, 'b');%-1.9
cs_model.modelNH4_Nlim_030 = modelNH4_Nlim_030;

%% NH4_Nlim_035
modelNH4_Nlim_035 = scaleBioMass(modelNH4_Nlim_035,'protein',0.49431);
modelNH4_Nlim_035 = scaleBioMass(modelNH4_Nlim_035,'RNA',0.07749);
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_2111',	0.35,	'b');
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1714',	-13.088,	'b');
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1654',	-1000,'l');
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1672',	exp1e(6,11),'b');%co2
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1992',	-exp1e(6,12),'b');%o2
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1761',	exp1e(6,13),'b');%ethanol
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1634',	exp1e(6,14),'b');%acetate 
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_2033',	exp1e(6,15),'b');%pyruvate 
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_2056',	exp1e(6,16),'b');%succinate 
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1808',	exp1e(6,17),'b');%glycerol 
modelNH4_Nlim_035 = changeObjective(modelNH4_Nlim_035,	'r_1654');
solutionNH4_Nlim_035 = optimizeCbModel(modelNH4_Nlim_035,'max');
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1654',	solutionNH4_Nlim_035.f, 'b');%-2.7
cs_model.modelNH4_Nlim_035 = modelNH4_Nlim_035;

%% C/N = 30
modelCN30 = scaleBioMass(modelCN30,'protein',0.3665);
modelCN30 = scaleBioMass(modelCN30,'RNA',0.04);
modelCN30 = scaleBioMass(modelCN30,'carbohydrate',0.4807);
modelCN30 = changeRxnBounds(modelCN30,	'r_2111',	0.2,	'b');
modelCN30 = changeRxnBounds(modelCN30,	'r_1714',	-1000,	'l');
modelCN30 = changeRxnBounds(modelCN30,	'r_1654',	-1000,	'l');
modelCN30 = changeRxnBounds(modelCN30,	'r_1672',	si4(2,3),'b');%co2
modelCN30 = changeRxnBounds(modelCN30,	'r_1992',	-si4(3,3),'b');%o2
modelCN30 = changeRxnBounds(modelCN30,	'r_1761',	si4(4,3),'b');%ethanol
modelCN30 = changeRxnBounds(modelCN30,	'r_1634',	si4(5,3),'b');%acetate 
modelCN30 = changeRxnBounds(modelCN30,	'r_2033',	si4(6,3),'b');%pyruvate 
modelCN30 = changeRxnBounds(modelCN30,	'r_2056',	si4(7,3),'b');%succinate 
modelCN30 = changeRxnBounds(modelCN30,	'r_1808',	si4(8,3),'b');%glycerol 
% GAM
modelCN30.S(324,3413) = -52.683;
modelCN30.S(610,3413) = -52.683;
modelCN30.S(288,3413) = 52.683;
modelCN30.S(601,3413) = 52.683;
modelCN30.S(1035,3413) = 52.683;
modelCN30 = changeObjective(modelCN30,	'r_1654');
solutionCN30 = optimizeCbModel(modelCN30,'max');
modelCN30 = changeRxnBounds(modelCN30,	'r_1654',	solutionCN30.f, 'b');%-1.157
cs_model.modelCN30 = modelCN30;

%% C/N = 50
modelCN50 = scaleBioMass(modelCN50,'protein',0.2635);
modelCN50 = scaleBioMass(modelCN50,'RNA',0.0195);
modelCN50 = scaleBioMass(modelCN50,'carbohydrate',0.6042);
modelCN50 = changeRxnBounds(modelCN50,	'r_2111',	0.2,	'b');
modelCN50 = changeRxnBounds(modelCN50,	'r_1714',	-1000,	'l');
modelCN50 = changeRxnBounds(modelCN50,	'r_1654',	-1000,	'l');
modelCN50 = changeRxnBounds(modelCN50,	'r_1672',	si4(2,4),'b');%co2
modelCN50 = changeRxnBounds(modelCN50,	'r_1992',	-si4(3,4),'b');%o2
modelCN50 = changeRxnBounds(modelCN50,	'r_1761',	si4(4,4),'b');%ethanol
modelCN50 = changeRxnBounds(modelCN50,	'r_1634',	si4(5,4),'b');%acetate 
modelCN50 = changeRxnBounds(modelCN50,	'r_2033',	si4(6,4),'b');%pyruvate 
modelCN50 = changeRxnBounds(modelCN50,	'r_2056',	si4(7,4),'b');%succinate 
modelCN50 = changeRxnBounds(modelCN50,	'r_1808',	si4(8,4),'b');%glycerol 
modelCN50.S(324,3413) = -49.902;
modelCN50.S(610,3413) = -49.902;
modelCN50.S(288,3413) = 49.902;
modelCN50.S(601,3413) = 49.902;
modelCN50.S(1035,3413) = 49.902;
modelCN50 = changeObjective(modelCN50,	'r_1654');
solutionCN50 = optimizeCbModel(modelCN50,'max');
modelCN50 = changeRxnBounds(modelCN50,	'r_1654',	solutionCN50.f, 'b');%-0.87
cs_model.modelCN50 = modelCN50;

%% C/N = 115
modelCN115 = scaleBioMass(modelCN115,'protein',0.2635);
modelCN115 = scaleBioMass(modelCN115,'RNA',0.0195);
modelCN115 = scaleBioMass(modelCN115,'carbohydrate',0.6042);
modelCN115 = changeRxnBounds(modelCN115,	'r_2111',	0.2,	'b');
modelCN115 = changeRxnBounds(modelCN115,	'r_1714',	-1000,	'l');
modelCN115 = changeRxnBounds(modelCN115,	'r_1654',	-1000,	'l');
modelCN115 = changeRxnBounds(modelCN115,	'r_1672',	si4(2,5),'b');%co2
modelCN115 = changeRxnBounds(modelCN115,	'r_1992',	-si4(3,5),'b');%o2
modelCN115 = changeRxnBounds(modelCN115,	'r_1761',	si4(4,5),'b');%ethanol
modelCN115 = changeRxnBounds(modelCN115,	'r_1634',	si4(5,5),'b');%acetate 
modelCN115 = changeRxnBounds(modelCN115,	'r_2033',	si4(6,5),'b');%pyruvate 
modelCN115 = changeRxnBounds(modelCN115,	'r_2056',	si4(7,5),'b');%succinate 
modelCN115 = changeRxnBounds(modelCN115,	'r_1808',	si4(8,5),'b');%glycerol 
modelCN115.S(324,3413) = -49.902;
modelCN115.S(610,3413) = -49.902;
modelCN115.S(288,3413) = 49.902;
modelCN115.S(601,3413) = 49.902;
modelCN115.S(1035,3413) = 49.902;
modelCN115 = changeObjective(modelCN115,	'r_1654');
solutionCN115 = optimizeCbModel(modelCN115,'max');
modelCN115 = changeRxnBounds(modelCN115,	'r_1654',	solutionCN115.f, 'b');%-0.86
cs_model.modelCN115 = modelCN115;

%% save model
cd output/Nlim_model
writeCbModel(cs_model.modelNH4_Nlim_005,'sbml','modelNH4_Nlim_005.xml')
writeCbModel(cs_model.modelNH4_Nlim_01,'sbml','modelNH4_Nlim_01.xml')
writeCbModel(cs_model.modelNH4_Nlim_013,'sbml','modelNH4_Nlim_013.xml')
writeCbModel(cs_model.modelNH4_Nlim_018,'sbml','modelNH4_Nlim_018.xml')
writeCbModel(cs_model.modelNH4_Nlim_030,'sbml','modelNH4_Nlim_030.xml')
writeCbModel(cs_model.modelNH4_Nlim_035,'sbml','modelNH4_Nlim_035.xml')
writeCbModel(cs_model.modelCN30,'sbml','modelCN30.xml')
writeCbModel(cs_model.modelCN50,'sbml','modelCN50.xml')
writeCbModel(cs_model.modelCN115,'sbml','modelCN115.xml')
writeCbModel(cs_model.modelgln_Nlim_01,'sbml','modelgln_Nlim_01.xml')
writeCbModel(cs_model.modelile_Nlim_01,'sbml','modelile_Nlim_01.xml')
writeCbModel(cs_model.modelphe_Nlim_01,'sbml','modelphe_Nlim_01.xml')
end