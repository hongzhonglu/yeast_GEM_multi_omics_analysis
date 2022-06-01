import pandas as pd
import cobra
import os
file = r'singlecell/output'
filename = os.listdir(file)
filepath = []
for f in filename:
    filepath.append(file + '/' + f)

model8 = cobra.io.read_sbml_model(r'singlecell/data/yeast-GEM.xml')
pfba_so = cobra.flux_analysis.pfba(model8)

rxnid = [i.id for i in model8.reactions]
expid = [e.split('.')[0] for e in filename]
flux_0 = [0.0]*len(rxnid)
pca_data = pd.DataFrame(index=rxnid, columns=expid)
flux = dict(zip(rxnid, flux_0))
c = 1
rxnnum = []
metnum = []

for s, en in zip(filepath, expid):
    print('processing model Number' + c)
    tempflux = flux.copy()
    tempmodel = cobra.io.read_sbml_model(s)
    temp_pfba = cobra.flux_analysis.pfba(tempmodel)
    temprxns = temp_pfba.fluxes.index.values.tolist()
    for tr in temprxns:
        tempflux[tr] = temp_pfba.fluxes.loc[tr]
    pca_data.loc[:, en] = tempflux.values()
    rxnnum.append(len(tempmodel.reactions))
    metnum.append(len(tempmodel.metabolites))

#for s, en in zip(filepath, expid):
#   print('processing model Number' + str(c))
#    tempmodel = cobra.io.read_sbml_model(s)
#    rxnnum.append(len(tempmodel.reactions))
#    metnum.append(len(tempmodel.metabolites))
#    c += 1
pca_data = pca_data.T
rrr = [r/20 for r in rxnnum]
mmm = [m/20 for m in metnum]
pca_data.insert(0, 'rxnnum/20', rrr)
pca_data.insert(0, 'metnum/20', mmm)
pca_data.to_excel('singlecell/output/pfba_flux.xlsx')

#pca_data = pca_data.drop('mmm', axis=1)
#pca_data = pca_data.drop('rrr', axis=1)

# 可视化
