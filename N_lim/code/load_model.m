% get N lim condition-specific model based on rosemarys' dataset
function [cs_model] = load_model
cd ..
model = importModel('data/yeast-GEM.xml');
% model = importModel('data/yeast_7.6_cobra.xml');
exp1e = readmatrix('data/1e mean.xlsx');
si4 = readmatrix('data/Supplementary Data 4.xlsx');
modelNH4_Nlim_01 = model; modelgln_Nlim_01 = model; modelphe_Nlim_01 = model; modelile_Nlim_01 = model;
modelNH4_Nlim_005 = model; modelNH4_Nlim_013 = model; modelNH4_Nlim_018 = model; modelNH4_Nlim_030 =model; modelNH4_Nlim_035 = model;
modelCN30 = model; modelCN50 = model; modelCN115 = model;
% scaleBioMass is a function in yeast-GEM
% yeast-GEM\code\otherChanges
cd ../../yeast-GEM/code/otherChanges/
%% NH4_Nlim_01
% constrain growth, exchange reaction, biomass
modelNH4_Nlim_01 = scaleBioMass(modelNH4_Nlim_01, 'protein', 0.22557);
modelNH4_Nlim_01 = scaleBioMass(modelNH4_Nlim_01, 'RNA', 0.0253);
% fix 90% measured growth rate as lower bound.
% 110% measured substrate uptake rate as lower bound.
% 90% measured outflow as lower bound.
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_2111',	0.1 * 0.8, 'l');
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1714',	-exp1e(2,10) * 1.2, 'l');%glucose
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1654',	-1000, 'l'); 
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1672',	exp1e(2,11) * 0.8,'l');%co2
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1992',	-exp1e(2,12) * 1.2 ,'l');%o2
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1761',	exp1e(2,13) * 0.8,'l');%ethanol
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1634',	exp1e(2,14) * 0.8,'l');%acetate 
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_2033',	exp1e(2,15) * 0.8,'l');%pyruvate 
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_2056',	exp1e(2,16) * 0.8,'l');%succinate 
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1808',	exp1e(2,17) * 0.8,'l');%glycerol 
modelNH4_Nlim_01 = changeObjective(modelNH4_Nlim_01,	'r_1654');
solutionNH4_Nlim_01 = optimizeCbModel(modelNH4_Nlim_01, 'max');
% constrain 110% minimized N usage
modelNH4_Nlim_01 = changeRxnBounds(modelNH4_Nlim_01,	'r_1654',	1.2 * solutionNH4_Nlim_01.f, 'l');
cs_model.modelNH4_Nlim_01 = modelNH4_Nlim_01;

%% gln_Nlim_01
modelgln_Nlim_01 = scaleBioMass(modelgln_Nlim_01,'protein',0.30677);
modelgln_Nlim_01 = scaleBioMass(modelgln_Nlim_01,'RNA',0.02215);
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_2111',	0.1 * 0.8,	'l');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1714',	-exp1e(9,10) * 1.2,	'l');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1654',	0,  'b');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1891',	-1000, 'l');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1672',	exp1e(9,11) * 0.8,'l');%co2
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1992',	-exp1e(9,12) * 1.2,'l');%o2
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1761',	exp1e(9,13) * 0.8,'l');%ethanol
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1634',	exp1e(9,14) * 0.8,'l');%acetate 
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_2033',	exp1e(9,15) * 0.8,'l');%pyruvate 
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_2056',	exp1e(9,16) * 0.8,'l');%succinate 
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1808',	exp1e(9,17) * 0.8,'l');%glycerol 
% modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1891',	-exp1e(9,9) * 1.2,'l');

modelgln_Nlim_01 = changeObjective(modelgln_Nlim_01,	'r_1891');
solutiongln_Nlim_01 = optimizeCbModel(modelgln_Nlim_01, 'max');
modelgln_Nlim_01 = changeRxnBounds(modelgln_Nlim_01,	'r_1891',	1.2 * solutiongln_Nlim_01.f, 'l');
cs_model.modelgln_Nlim_01 = modelgln_Nlim_01;

%% phe_Nlim_01
modelphe_Nlim_01 = scaleBioMass(modelphe_Nlim_01,'protein',0.3985);
modelphe_Nlim_01 = scaleBioMass(modelphe_Nlim_01,'RNA',0.02423);
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_2111',	0.1 * 0.8,	'l');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1714',	-exp1e(11,10) * 1.2, 'l');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1654',	0,'b');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1903',	-1000,'l');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1672',	exp1e(11,11) * 0.8,'l');%co2
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1992',	-exp1e(11,12) * 1.2,'l');%o2
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1761',	exp1e(11,13) * 0.8,'l');%ethanol
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1634',	exp1e(11,14) * 0.8,'l');%acetate 
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_2033',	exp1e(11,15) * 0.8,'l');%pyruvate 
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_2056',	exp1e(11,16) * 0.8,'l');%succinate 
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1808',	exp1e(11,17) * 0.8,'l');%glycerol 
%modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1903',	-exp1e(11,9) * 1.2, 'l');

modelphe_Nlim_01 = changeObjective(modelphe_Nlim_01,	'r_1903');
solutionphe_Nlim_01 = optimizeCbModel(modelphe_Nlim_01, 'max');
modelphe_Nlim_01 = changeRxnBounds(modelphe_Nlim_01,	'r_1903',	solutionphe_Nlim_01.f * 1.2, 'l');
cs_model.modelphe_Nlim_01 = modelphe_Nlim_01;

%% ile_Nlim_01
modelile_Nlim_01 = scaleBioMass(modelile_Nlim_01,'protein',0.43732);
modelile_Nlim_01 = scaleBioMass(modelile_Nlim_01,'RNA',0.02381);
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_2111',	0.1 * 0.8,	'l');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1714',	-exp1e(13,10) * 1.2,	'l');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1654',	0,'b');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1897',	-1000,'l');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1672',	exp1e(13,11) * 0.8,'l');%co2
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1992',	-exp1e(13,12) * 1.2,'l');%o2
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1761',	exp1e(13,13) * 0.8,'l');%ethanol
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1634',	exp1e(13,14) * 0.8,'l');%acetate 
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_2033',	exp1e(13,15) * 0.8,'l');%pyruvate 
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_2056',	exp1e(13,16) * 0.8,'l');%succinate
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1808',	exp1e(13,17) * 0.8,'l');%glycerol 
modelile_Nlim_01 = changeObjective(modelile_Nlim_01,	'r_1897');
solutionile_Nlim_01 = optimizeCbModel(modelile_Nlim_01, 'max');
modelile_Nlim_01 = changeRxnBounds(modelile_Nlim_01,	'r_1897',	solutionile_Nlim_01.f * 1.2, 'l');
cs_model.modelile_Nlim_01 = modelile_Nlim_01;

%% NH4_Nlim_005
modelNH4_Nlim_005 = scaleBioMass(modelNH4_Nlim_005,'protein',0.26148);
modelNH4_Nlim_005 = scaleBioMass(modelNH4_Nlim_005,'RNA',0.02417);
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_2111',	0.05 * 0.8,	'l');
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1714',	-exp1e(1,10) * 1.2,	'l');
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1654',	-1000,'l');
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1672',	exp1e(1,11) * 0.8,'l');%co2
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1992',	-exp1e(1,12) * 1.2,'l');%o2
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1761',	exp1e(1,13) * 0.8,'l');%ethanol
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1634',	exp1e(1,14) * 0.8,'l');%acetate 
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_2033',	exp1e(1,15) * 0.8,'l');%pyruvate 
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_2056',	exp1e(1,16) * 0.8,'l');%succinate 
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1808',	exp1e(1,17) * 0.8,'l');%glycerol 
modelNH4_Nlim_005 = changeObjective(modelNH4_Nlim_005,	'r_1654');
solutionNH4_Nlim_005 = optimizeCbModel(modelNH4_Nlim_005, 'max');
modelNH4_Nlim_005 = changeRxnBounds(modelNH4_Nlim_005,	'r_1654',	solutionNH4_Nlim_005.f * 1.2, 'l');
cs_model.modelNH4_Nlim_005 = modelNH4_Nlim_005;

%% NH4_Nlim_013
modelNH4_Nlim_013 = scaleBioMass(modelNH4_Nlim_013,'protein',0.27569);
modelNH4_Nlim_013 = scaleBioMass(modelNH4_Nlim_013,'RNA',0.02583);
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_2111',	0.13 * 0.8,	'l');
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1714',	-exp1e(3,10) * 1.2,	'l');
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1654',	-1000,'l');
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1672',	exp1e(3,11) * 0.8,'l');%co2
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1992',	-exp1e(3,12) * 1.2,'l');%o2
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1761',	exp1e(3,13) * 0.8,'l');%ethanol
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1634',	exp1e(3,14) * 0.8,'l');%acetate 
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_2033',	exp1e(3,15) * 0.8,'l');%pyruvate 
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_2056',	exp1e(3,16) * 0.8,'l');%succinate 
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1808',	exp1e(3,17) * 0.8,'l');%glycerol 
modelNH4_Nlim_013 = changeObjective(modelNH4_Nlim_013,	'r_1654');
solutionNH4_Nlim_013 = optimizeCbModel(modelNH4_Nlim_013, 'max');
modelNH4_Nlim_013 = changeRxnBounds(modelNH4_Nlim_013,	'r_1654',	solutionNH4_Nlim_013.f * 1.2, 'l');
cs_model.modelNH4_Nlim_013 = modelNH4_Nlim_013;

%% NH4_Nlim_018
modelNH4_Nlim_018 = scaleBioMass(modelNH4_Nlim_018,'protein',0.28647);
modelNH4_Nlim_018 = scaleBioMass(modelNH4_Nlim_018,'RNA',0.03854);
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_2111',	0.18 * 0.8,	'l');
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1714',	-exp1e(4,10) * 1.2, 'l');
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1654',	-1000,'l');
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1672',	exp1e(4,11) * 0.8,'l');%co2
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1992',	-exp1e(4,12) * 1.2,'l');%o2
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1761',	exp1e(4,13) * 0.8,'l');%ethanol
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1634',	exp1e(4,14) * 0.8,'l');%acetate 
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_2033',	exp1e(4,15) * 0.8,'l');%pyruvate 
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_2056',	exp1e(4,16) * 0.8,'l');%succinate 
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1808',	exp1e(4,17) * 0.8,'l');%glycerol 
modelNH4_Nlim_018 = changeObjective(modelNH4_Nlim_018,	'r_1654');
solutionNH4_Nlim_018 = optimizeCbModel(modelNH4_Nlim_018, 'max');
modelNH4_Nlim_018 = changeRxnBounds(modelNH4_Nlim_018,	'r_1654',	solutionNH4_Nlim_018.f * 1.2, 'l');%-0.86
cs_model.modelNH4_Nlim_018 = modelNH4_Nlim_018;

%% NH4_Nlim_030
modelNH4_Nlim_030 = scaleBioMass(modelNH4_Nlim_030,'protein',0.41448);
modelNH4_Nlim_030 = scaleBioMass(modelNH4_Nlim_030,'RNA',0.04822);
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_2111',	0.30 * 0.8,	'l');
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1714',	-exp1e(5,10) * 1.2,	'l');
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1654',	-1000,'l');
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1672',	exp1e(5,11) * 0.8,'l');%co2
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1992',	-exp1e(5,12) * 1.2,'l');%o2
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1761',	exp1e(5,13) * 0.8,'l');%ethanol
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1634',	exp1e(5,14) * 0.8,'l');%acetate 
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_2033',	exp1e(5,15) * 0.8,'l');%pyruvate 
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_2056',	exp1e(5,16) * 0.8,'l');%succinate 
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1808',	exp1e(5,17) * 0.8,'l');%glycerol 
modelNH4_Nlim_030 = changeObjective(modelNH4_Nlim_030,	'r_1654');
solutionNH4_Nlim_030 = optimizeCbModel(modelNH4_Nlim_030, 'max');
modelNH4_Nlim_030 = changeRxnBounds(modelNH4_Nlim_030,	'r_1654',	solutionNH4_Nlim_030.f * 1.2, 'l');
cs_model.modelNH4_Nlim_030 = modelNH4_Nlim_030;

%% NH4_Nlim_035
modelNH4_Nlim_035 = scaleBioMass(modelNH4_Nlim_035,'protein',0.49431);
modelNH4_Nlim_035 = scaleBioMass(modelNH4_Nlim_035,'RNA',0.07749);
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_2111',	0.35 * 0.8,	'l');
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1714',	-exp1e(6,10) * 1.2,	'l');
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1654',	-1000,'l');
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1672',	exp1e(6,11) * 0.8,'l');%co2
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1992',	-exp1e(6,12) * 1.2,'l');%o2
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1761',	exp1e(6,13) * 0.8,'l');%ethanol
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1634',	exp1e(6,14) * 0.8,'l');%acetate 
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_2033',	exp1e(6,15) * 0.8,'l');%pyruvate 
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_2056',	exp1e(6,16) * 0.8,'l');%succinate 
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1808',	exp1e(6,17) * 0.8,'l');%glycerol 
modelNH4_Nlim_035 = changeObjective(modelNH4_Nlim_035,	'r_1654');
solutionNH4_Nlim_035 = optimizeCbModel(modelNH4_Nlim_035, 'max');
modelNH4_Nlim_035 = changeRxnBounds(modelNH4_Nlim_035,	'r_1654',	solutionNH4_Nlim_035.f * 1.2, 'l');
cs_model.modelNH4_Nlim_035 = modelNH4_Nlim_035;

%% C/N = 30
modelCN30 = scaleBioMass(modelCN30,'protein',0.3665);
modelCN30 = scaleBioMass(modelCN30,'RNA',0.04);
modelCN30 = scaleBioMass(modelCN30,'carbohydrate',0.4807);
modelCN30 = changeRxnBounds(modelCN30,	'r_2111',	0.2 * 0.8,	'l');
modelCN30 = changeRxnBounds(modelCN30,	'r_1714',	-si4(1,3) * 1.2,	'l');
modelCN30 = changeRxnBounds(modelCN30,	'r_1654',	-1000,	'l');
modelCN30 = changeRxnBounds(modelCN30,	'r_1672',	si4(2,3) * 0.8,'l');%co2
modelCN30 = changeRxnBounds(modelCN30,	'r_1992',	-si4(3,3) * 1.2,'l');%o2
modelCN30 = changeRxnBounds(modelCN30,	'r_1761',	si4(4,3) * 0.8,'l');%ethanol
modelCN30 = changeRxnBounds(modelCN30,	'r_1634',	si4(5,3) * 0.8,'l');%acetate 
modelCN30 = changeRxnBounds(modelCN30,	'r_2033',	si4(6,3) * 0.8,'l');%pyruvate 
modelCN30 = changeRxnBounds(modelCN30,	'r_2056',	si4(7,3) * 0.8,'l');%succinate 
modelCN30 = changeRxnBounds(modelCN30,	'r_1808',	si4(8,3) * 0.8,'l');%glycerol 
% GAM
%modelCN30.S(324,3413) = -52.683;
%modelCN30.S(610,3413) = -52.683;
%modelCN30.S(288,3413) = 52.683;
%modelCN30.S(601,3413) = 52.683;
%modelCN30.S(1035,3413) = 52.683;
modelCN30 = changeObjective(modelCN30,	'r_1654');
solutionCN30 = optimizeCbModel(modelCN30, 'max');
modelCN30 = changeRxnBounds(modelCN30,	'r_1654',	1.2 * solutionCN30.f, 'l');
cs_model.modelCN30 = modelCN30;

%% C/N = 50
modelCN50 = scaleBioMass(modelCN50,'protein',0.2635);
modelCN50 = scaleBioMass(modelCN50,'RNA',0.0195);
modelCN50 = scaleBioMass(modelCN50,'carbohydrate',0.6042);
modelCN50 = changeRxnBounds(modelCN50,	'r_2111',	0.2 * 0.8,	'l');
modelCN50 = changeRxnBounds(modelCN50,	'r_1714',	-si4(1,4) * 1.2,	'l');
modelCN50 = changeRxnBounds(modelCN50,	'r_1654',	-1000,	'l');
modelCN50 = changeRxnBounds(modelCN50,	'r_1672',	si4(2,4) * 0.8,'l');%co2
modelCN50 = changeRxnBounds(modelCN50,	'r_1992',	-si4(3,4) * 1.2,'l');%o2
modelCN50 = changeRxnBounds(modelCN50,	'r_1761',	si4(4,4) * 0.8,'l');%ethanol
modelCN50 = changeRxnBounds(modelCN50,	'r_1634',	si4(5,4) * 0.8,'l');%acetate 
modelCN50 = changeRxnBounds(modelCN50,	'r_2033',	si4(6,4) * 0.8,'l');%pyruvate 
modelCN50 = changeRxnBounds(modelCN50,	'r_2056',	si4(7,4) * 0.8,'l');%succinate 
modelCN50 = changeRxnBounds(modelCN50,	'r_1808',	si4(8,4) * 0.8,'l');%glycerol 
%modelCN50.S(324,3413) = -49.902;
%modelCN50.S(610,3413) = -49.902;
%modelCN50.S(288,3413) = 49.902;
%modelCN50.S(601,3413) = 49.902;
%modelCN50.S(1035,3413) = 49.902;
modelCN50 = changeObjective(modelCN50,	'r_1654');
solutionCN50 = optimizeCbModel(modelCN50, 'max');
modelCN50 = changeRxnBounds(modelCN50,	'r_1654',	1.2 * solutionCN50.f, 'l');
cs_model.modelCN50 = modelCN50;

%% C/N = 115
modelCN115 = scaleBioMass(modelCN115,'protein',0.2635);
modelCN115 = scaleBioMass(modelCN115,'RNA',0.0195);
modelCN115 = scaleBioMass(modelCN115,'carbohydrate',0.6042);
modelCN115 = changeRxnBounds(modelCN115,	'r_2111',	0.2 * 0.8,	'l');
modelCN115 = changeRxnBounds(modelCN115,	'r_1714',	-si4(1,5) * 1.2,	'l');
modelCN115 = changeRxnBounds(modelCN115,	'r_1654',	-1000,	'l');
modelCN115 = changeRxnBounds(modelCN115,	'r_1672',	si4(2,5) * 0.8,'l');%co2
modelCN115 = changeRxnBounds(modelCN115,	'r_1992',	-si4(3,5) * 1.2,'l');%o2
modelCN115 = changeRxnBounds(modelCN115,	'r_1761',	si4(4,5) * 0.8,'l');%ethanol
modelCN115 = changeRxnBounds(modelCN115,	'r_1634',	si4(5,5) * 0.8,'l');%acetate 
modelCN115 = changeRxnBounds(modelCN115,	'r_2033',	si4(6,5) * 0.8,'l');%pyruvate 
modelCN115 = changeRxnBounds(modelCN115,	'r_2056',	si4(7,5) * 0.8,'l');%succinate 
modelCN115 = changeRxnBounds(modelCN115,	'r_1808',	si4(8,5) * 0.8,'l');%glycerol 
%modelCN115.S(324,3413) = -49.902;
%modelCN115.S(610,3413) = -49.902;
%modelCN115.S(288,3413) = 49.902;
%modelCN115.S(601,3413) = 49.902;
%modelCN115.S(1035,3413) = 49.902;
modelCN115 = changeObjective(modelCN115,	'r_1654');
solutionCN115 = optimizeCbModel(modelCN115, 'max');
modelCN115 = changeRxnBounds(modelCN115,	'r_1654',	1.2 * solutionCN115.f, 'l');
cs_model.modelCN115 = modelCN115;

%% save model
cd ../../../yeast_GEM_multi_omics_analysis/N_lim/output/Nlim_model/
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