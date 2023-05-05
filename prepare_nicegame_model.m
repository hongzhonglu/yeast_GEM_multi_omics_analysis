seeddata = readtable('D:\model_research\yeast_GEM_multi_omics_analysis\seeddata.xlsx');
yeast = importModel("D:\model_research\yeast-GEM\model\yeast-GEM.xml");
% newmodel.metCompSymbol = model.metCompSymbol
yeast.metCompSymbol = cell(length(yeast.metNames), 1);
for i = 1:length(yeast.metComps)
    yeast.metCompSymbol(i, 1) = yeast.comps(yeast.metComps(i));
end
% newmodel.metSEEDID = model.metSEEDID
na = 0;
for m = 1:length(yeast.metNames)

    try
        met_kegg = yeast.metMiriams{m}.value(find(strcmp('kegg.compound', yeast.metMiriams{m}.name)));
        met_kegg = met_kegg{1};
    catch
        met_kegg = 'NA';
        na = na + 1;
    end

    yeast.metSEEDID{m, 1} = met_kegg;
end
% newmodel.CompartmentData = model.CompartmentData
yeast.CompartmentData.compSymbolList = yeast.comps.';
yeast.CompartmentData.compNameList = yeast.compNames.';
yeast.CompartmentData.pH = [7.5,7,7,7,7,7,7,7,7,7,7,7,7,7];
yeast.CompartmentData.ioninStr = [0.25,0,0,0,0,0,0,0,0,0,0,0,0,0];
yeast.CompartmentData.compMaxConc = [0.0500,0.0500,0.1000,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500,0.0500];
yeast.CompartmentData.compMinConc = [1.000e-06,1.000e-06,1.000e-08,1.000e-06,1.000e-06,1.000e-06,1.000e-06,1.000e-06,1.000e-06,1.000e-06,1.000e-06,1.000e-06,1.000e-06,1.000e-06];
yeast.CompartmentData.membranePot = zeros(14, 14);
yeast.CompartmentData.membranePot(14, 1) = -150;
yeast.CompartmentData.membranePot(1, 14) = 150;
