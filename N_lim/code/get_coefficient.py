###########################################################################
# Description
# get coefficient of the chosen reaction
# ignore 0 flux and protein more than 7 out of 12 condition of the reaction
###########################################################################
# Parameter
# flux_prot: A dictionary. The output of 'paring'.
###########################################################################
# Output
# slope: coefficient of the chosen reaction.
###########################################################################

import numpy as np
import pandas as pd
import scipy.stats as st
def get_coefficient(flux_prot):
    conditionname = ['CN115', 'CN50', 'CN30', 'Phe', 'Ile',
                     'N035', 'N030', 'N018', 'N013',
                     'N010', 'N005', 'Gln']
    # get positive flux
    fluxtemp = flux_prot.loc[:, 'flux'].values.tolist()
    fluxtemp2 = [abs(f) for f in fluxtemp]
    # filtrate 0 flux
    fluxvalue = [fl if fl > 0.000001 else 0 for fl in fluxtemp2]
    prottemp = flux_prot.loc[:, 'prot'].values.tolist()
    if len(str(prottemp[0]).split(' or ')) == 1:
        # no 'or'
        protvalue = prottemp
        slope = linear_regression(fluxvalue, protvalue)
    else:
        # isoenzyme
        split_pro = [s.split(' or ') for s in prottemp]
        slopelist = []
        for l in range(len(split_pro[0])):
            protvalue = []
            for sp in split_pro:
                protvalue.append(float(sp[l]))
            tempslope = linear_regression(fluxvalue, protvalue)
            if tempslope != 'nan':
                slopelist.append(tempslope)
        if len(slopelist) != 0:
            slope = np.mean(slopelist)
        else:
            slope = 'nan'
    return slope


def get_coefficient_z(flux_prot):
    conditionname = ['CN115', 'CN50', 'CN30', 'Phe', 'Ile',
                     'N035', 'N030', 'N018', 'N013',
                     'N010', 'N005', 'Gln']
    # get positive flux
    fluxvalue = flux_prot.loc[:, 'flux'].values.tolist()
    prottemp = flux_prot.loc[:, 'prot'].values.tolist()
    if len(str(prottemp[0]).split(' or ')) == 1:
        # no 'or'
        protvalue = prottemp
        slope = pearson_zscore(fluxvalue, protvalue)
    else:
        # isoenzyme
        split_pro = [s.split(' or ') for s in prottemp]
        slopelist = []
        for l in range(len(split_pro[0])):
            protvalue = []
            for sp in split_pro:
                protvalue.append(float(sp[l]))
            tempslope = pearson_zscore(fluxvalue, protvalue)
            if tempslope != 'nan':
                slopelist.append(tempslope)
        if len(slopelist) != 0:
            slope = np.mean(slopelist)
        else:
            slope = 'nan'
    return slope


# linear regression to get slope
def linear_regression(fluxvalue, protvalue):
    fluxvalue = [float(ff) for ff in fluxvalue]
    protvalue = [float(pp) for pp in protvalue]
    nozerofluxindex = np.nonzero(fluxvalue)
    nozeroprotindex = np.nonzero(protvalue)
    commonnonzero = np.intersect1d(nozerofluxindex, nozeroprotindex)
    if len(commonnonzero) > 5:
        flux = [fluxvalue[flu] for flu in commonnonzero]
        protein = [protvalue[pro] for pro in commonnonzero]
        logflux = np.log10(flux)
        logprot = np.log10(protein)
        slope, intercept, r_value, p_value, std_err = st.linregress(logflux, logprot)
    else:
        slope = 'nan'
    return slope

def pearson_zscore(fluxvalue, protvalue):
    fluxvalue = [float(ff) for ff in fluxvalue]
    protvalue = [float(pp) for pp in protvalue]
    # data = pd.DataFrame(columns=['f', 'p'])
    # data.loc[:, 'f'] = fluxvalue
    # data.loc[:, 'f'] = fluxvalue
    if len(set(fluxvalue)) > 1 or len(set(protvalue)) > 1:
        rho = np.corrcoef(fluxvalue, protvalue)
    else:
        rho = 'nan'
    return rho[1, 0]
