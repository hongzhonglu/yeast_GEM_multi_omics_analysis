###########################################################################
# Description
# paring reactionID with its corresponding flux and enzyme concentration
# for isoenzyme, select all the concentration.
# for complex, select the lowest protein concentration as the complex concentration
###########################################################################
# Parameter
# fluxdict: A dictionary.
#           Keys are condition name,
#           Values are dataframe containing reactionID(index), flux and gene name
# protdict: A dictionary.
#           Keys are condition name,
#           Values are dataframe containing gene name(index) and protein concentration
# reactionID: Assigned reaction by user
###########################################################################
# Output
# flux_prot: A dataframe.
#            index is condition name
#            columns is flux and protein concentration of the corresponding reactionID
###########################################################################
import pandas as pd
import re
import numpy as np

def get_prot_conc(tempgene, sgd, protdict, condition):
    gene1 = re.compile(r'Y.{6}').findall(tempgene)
    gene2 = re.compile(r'Y.{6}-\w').findall(tempgene)
    gene3 = re.compile(r'Q.{4}').findall(tempgene)
    gene = gene1 + gene2 + gene3
    templist = []
    for t in gene:
        try:
            t = sgd.loc[t, 'Standard_name']
            templist.append(protdict[condition].loc[t, 'prot'])
        except KeyError or ValueError:
            pass
    return templist


def oneone(templist, flux_prot, condition):
    if len(templist) != 0:
        flux_prot.loc[condition, 'prot'] = templist[0]
    else:
        flux_prot.loc[condition, 'prot'] = 0
    return flux_prot


def complex(templist, flux_prot, condition):
    if len(templist) != 0:
        conc = np.min(templist)
        flux_prot.loc[condition, 'prot'] = np.min(templist)
    else:
        conc = [0]
        flux_prot.loc[condition, 'prot'] = 0
    return flux_prot, conc


def isoenzyme(templist, flux_prot, condition):
    flux_prot.loc[condition, 'prot'] = np.sum(templist)
    return flux_prot, np.sum(templist)


def paring(fluxdict, protdict, reactionID):
    sgd = pd.read_table(r'D:\model_research\yeast_GEM_multi_omics_analysis\N_lim\data\SGDgeneNames.tsv')
    sgd.index = sgd.loc[:, 'Systematic_name']
    conditionname = ['CN115', 'CN50', 'CN30', 'Phe', 'Ile',
                     'N035', 'N030', 'N018', 'N013',
                     'N010', 'N005', 'Gln']
    flux_prot = pd.DataFrame(index=conditionname,
                             columns=['flux', 'prot'])
    print(reactionID)
    for rr, cc in zip(fluxdict.keys(), conditionname):
        tempgene = fluxdict[rr].loc[reactionID, 'gene']
        re_and = re.compile('and')
        re_or = re.compile('or')
        com_and = re_and.findall(tempgene)
        com_or = re_or.findall(tempgene)
        flux_prot.loc[cc, 'flux'] = fluxdict[rr].loc[reactionID, 'flux']
        # one enzyme catalyse one reaction
        if len(com_or) == 0 and len(com_and) == 0:
            templist = get_prot_conc(tempgene, sgd, protdict, cc)
            flux_prot = oneone(templist, flux_prot, cc)
        # complex
        if len(com_or) == 0 and len(com_and) != 0:
            templist = get_prot_conc(tempgene, sgd, protdict, cc)
            flux_prot, __ = complex(templist, flux_prot, cc)
        # isoenzymes
        if len(com_or) != 0 and len(com_and) == 0:
            templist = get_prot_conc(tempgene, sgd, protdict, cc)
            flux_prot, __ = isoenzyme(templist, flux_prot, cc)
        # complex and isoenzyme
        if len(com_or) != 0 and len(com_and) != 0:
            tempcomplex = tempgene.split(' or ')
            isocom = []
            for tc in tempcomplex:
                templist = get_prot_conc(tc, sgd, protdict, cc)
                __, conc = complex(templist, flux_prot, cc)
                isocom = isocom + conc
            flux_prot.loc[cc, 'prot'] = np.sum(isocom)
    return flux_prot
