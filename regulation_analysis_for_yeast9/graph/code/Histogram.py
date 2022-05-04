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

def histogram(reaction_slope, savename):
    data = [d for d in reaction_slope.values() if d != 'nan']
    sns.set(style="dark")
    sns.distplot(data,
                 kde=True,
                 hist=True,
                 norm_hist=True,
                 kde_kws={"color": "g", "alpha": 0.3, "linewidth": 5, "shade": True})
    # in the next version of the distplot function, one would have to write:
    # sns.distplot(data=df, x="sepal_length", kde=True, kde_kws={"color": "g", "alpha": 0.3, "linewidth": 5, "shade": True})
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 16,
             }
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 20,
             }
    plt.xlabel('ρ', font1)
    plt.xlabel('Density', font1)
    plt.title(savename, font2)
    s_name = savename + '.png'
    plt.savefig(s_name)
    plt.show()
