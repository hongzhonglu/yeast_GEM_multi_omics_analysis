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

def paring(fluxdict, protdict, reactionID):
    conditionname = ['CN115', 'CN50', 'CN30', 'Phe', 'Ile',
                     'N035', 'N030', 'N018', 'N013',
                     'N010', 'N005', 'Gln']
    flux_prot = pd.DataFrame(index=conditionname,
                             columns=['flux', 'prot']
                            )
    N_protgene_list = protdict['Phe'].index.values.tolist()
    CN_protgene_list = protdict['CN30'].index.values.tolist()
    commongene = [val for val in N_protgene_list if val in CN_protgene_list]
    num = -1
    for rr, cc in zip(fluxdict.keys(), conditionname):
        num += 1
        tempgene = fluxdict[rr].loc[reactionID, 'gene']
        re_and = re.compile('and')
        re_or = re.compile('or')
        com_and = re_and.findall(tempgene)
        com_or = re_or.findall(tempgene)
        flux_prot.loc[cc, 'flux'] = fluxdict[rr].loc[reactionID, 'flux']
        # one enzyme catalyse one reaction
        if len(com_or) == 0 and len(com_and) == 0:
            if num < 3:
                try:
                    tempindex = commongene.index(tempgene)
                    flux_prot.loc[cc, 'prot'] = protdict[cc].loc[tempgene, 'prot']
                except ValueError or KeyError:
                    flux_prot.loc[cc, 'prot'] = 0
            else:
                try:
                    tempindex = commongene.index(tempgene)
                    flux_prot.loc[cc, 'prot'] = protdict[cc].loc[tempgene, 'prot']
                except ValueError or KeyError:
                    flux_prot.loc[cc, 'prot'] = 0
        # complex
        elif len(com_or) == 0 and len(com_and) != 0:
            if num < 3:
                tempgene = tempgene.split(' and ')
                templist = []
                for t in tempgene:
                    try:
                        tempindex = commongene.index(t)
                        templist.append(protdict[cc].loc[t, 'prot'])
                    except ValueError or KeyError:
                        templist.append(0)
                flux_prot.loc[cc, 'prot'] = min(templist)
            else:
                tempgene = tempgene.split(' and ')
                templist = []
                for t in tempgene:
                    try:
                        tempindex = commongene.index(t)
                        templist.append(protdict[cc].loc[t, 'prot'])
                    except ValueError or KeyError:
                        templist.append(0)
                flux_prot.loc[cc, 'prot'] = min(templist)
        # isoenzymes
        elif len(com_or) != 0 and len(com_and) == 0:
            if num < 3:
                tempgene = tempgene.split(' or ')
                templist = []
                for t in tempgene:
                    try:
                        tempindex = commongene.index(t)
                        templist.append(str(protdict[cc].loc[t, 'prot']))
                    except ValueError or KeyError:
                        templist.append('0')
                tempstr = ' or '.join(templist)
                flux_prot.loc[cc, 'prot'] = tempstr
            else:
                tempgene = tempgene.split(' or ')
                templist = []
                for t in tempgene:
                    try:
                        tempindex = commongene.index(t)
                        templist.append(str(protdict[cc].loc[t, 'prot']))
                    except ValueError or KeyError:
                        templist.append('0')
                tempstr = ' or '.join(templist)
                flux_prot.loc[cc, 'prot'] = tempstr

        # complex and isoenzyme
        elif len(com_or) != 0 and len(com_and) != 0:
            if num < 3:
                tempcomplex = tempgene.split(' or ')
                templist = []
                tempcomplexgene = []

                for t in tempcomplex:
                    if len(t.split('(')) == 1:
                        tempcomplexgene.append(t)
                    else:
                        t = t.split('(')[1].split(')')[0]
                        tempcomplexgene.append(t)
                if len(re_and.findall(tempcomplexgene[0])) == 0:
                    tempcomplexgene = tempcomplexgene[::-1]
                for tt in tempcomplexgene:
                    co_and = re_and.findall(tt)
                    if len(co_and) == 0:
                        try:
                            tempindex = commongene.index(tt)
                            templist.append(protdict[cc].loc[tt, 'prot'])
                        except ValueError or KeyError:
                            templist.append(0)
                    else:
                        tempcomlist = []
                        tt = tt.split(' and ')
                        for omg in tt:
                            try:
                                tempindex = commongene.index(omg)
                                tempcomlist.append(protdict[cc].loc[omg, 'prot'])
                            except ValueError or KeyError:
                                tempcomlist.append(0)
                    templist.append(min(tempcomlist))
                for te in range(len(templist)):
                    templist[te] = str(templist[te])
                tempstr = ' or '.join(templist)
                flux_prot.loc[cc, 'prot'] = tempstr
            else:
                tempcomplex = tempgene.split(' or ')
                templist = []
                tempcomplexgene = []

                for t in tempcomplex:
                    if len(t.split('(')) == 1:
                        tempcomplexgene.append(t)
                    else:
                        t = t.split('(')[1].split(')')[0]
                        tempcomplexgene.append(t)
                if len(re_and.findall(tempcomplexgene[0])) == 0:
                    tempcomplexgene = tempcomplexgene[::-1]
                for tt in tempcomplexgene:
                    co_and = re_and.findall(tt)
                    if len(co_and) == 0:
                        try:
                            tempindex = commongene.index(tt)
                            templist.append(protdict[cc].loc[tt, 'prot'])
                        except ValueError or KeyError:
                            templist.append(0)
                    else:
                        tempcomlist = []
                        tt = tt.split(' and ')
                        for omg in tt:
                            try:
                                tempindex = commongene.index(omg)
                                tempcomlist.append(protdict[cc].loc[omg, 'prot'])
                            except ValueError or KeyError:
                                tempcomlist.append(0)
                    templist.append(min(tempcomlist))
                for te in range(len(templist)):
                    templist[te] = str(templist[te])
                tempstr = ' or '.join(templist)
                flux_prot.loc[cc, 'prot'] = tempstr
    return flux_prot
