###########################################################################
# Description
# get histogram
###########################################################################
# Parameter
# reaction_slope: A dictionary.
#                 Keys are reactionIDs,
#                 Values are slopes.
# savename:       file name you wanna save.
###########################################################################
# Output
# histogram figure
###########################################################################

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def histogram(reaction_slope):
    lable = ['All reactions',
             'Carbohydrate metabolism',
             'Amino acid metbaolism',
             'Lipid metabolism']
    plt.figure(figsize=(6, 4))
    for i in [1, 2, 3, 4]:
        plt.subplot(2, 2, i)
        data = [reaction_slope[i-1].loc[d, 'p'] for d in reaction_slope[i-1].index if reaction_slope[i-1].loc[d, 'p'] != 'nan']
        colors = ['#FA499A', '#DE40D7', '#CD54F5', '#8C40DE']
        sns.distplot(data,
                     color='#C9ADF3',
                     kde_kws={"color": colors[i-1], "alpha": 0.2, "linewidth": 1, "shade": True})
        font1 = {'family': 'Arial',
                 'weight': 'normal',
                 'size': 12,
                 }
        plt.rc('font', **font1)
        plt.xlabel(lable[i-1], fontdict=font1)
        plt.xlim(xmax=1.8, xmin=-1.8)
        plt.ylabel('Density', fontdict=font1)
        if i == 1 or i == 2:
            plt.xticks([])
        if i == 2 or i == 4:
            plt.yticks([])
            plt.ylabel('')
        plt.tight_layout()
    plt.savefig('./output/p distri.tif', dpi=600)
    plt.show()
