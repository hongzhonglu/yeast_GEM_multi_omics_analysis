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

import matplotlib.pyplot as plt

from matplotlib.pyplot import MultipleLocator


def histogram(reaction_slope):
    lable = ['All reactions',
             'Carbohydrate metabolism',
             'Amino acid metbaolism',]

    plt.figure(figsize=(6, 3))
    for i in [1, 2, 3]:
        plt.subplot(1, 3, i)
        data = [reaction_slope[i-1].loc[d, 'p'] for d in reaction_slope[i-1].index if reaction_slope[i-1].loc[d, 'p'] != 'nan']
        colors = ['#FA499A', '#DE40D7', '#CD54F5', '#8C40DE']
        # sns.distplot(data,
        #              color='#C9ADF3',
        #              kde_kws={"color": colors[i-1], "alpha": 0.2, "linewidth": 1, "shade": True})
        if i == 1:
            plt.hist(data, bins=10, stacked=True, color='#f8cb7f')
        if i == 2:
            plt.hist(data, bins=6, stacked=True, color='#7cd6cf')
        if i == 3:
            plt.hist(data, bins=10, stacked=True, color='#7898e1')
            y_major_locator = MultipleLocator(3)
            ax = plt.gca()
            ax.yaxis.set_major_locator(y_major_locator)
        # if i == 4:
        #     plt.hist(data, bins=5, stacked=True, color='#9987ce')
        # if i == 1:
        #     plt.hist(data, bins=20, stacked=True, color='#f8cb7f')
        # if i == 2:
        #     plt.hist(data, bins=10, stacked=True, color='#7cd6cf')
        # if i == 3:
        #     plt.hist(data, bins=15, stacked=True, color='#7898e1')
        # if i == 4:
        #     plt.hist(data, bins=5, stacked=True, color='#9987ce')

        font1 = {'family': 'Arial',
                 'weight': 'normal',
                 'size': 18,
                 }
        plt.rc('font', **font1)
        plt.xlabel(lable[i-1], fontdict=font1)
        plt.xlim(xmax=1.8, xmin=-1.8)
        plt.ylabel('Reaction number', fontdict=font1)
        plt.yticks(fontproperties='Arial', size=18)
        plt.xlabel('ρ value', fontdict=font1)
        # if i == 1 or i == 2:
        #     plt.xticks([])
        # if i == 2 or i == 4:
        #     # plt.yticks([])
        #     plt.ylabel('')
    plt.tight_layout()
    plt.savefig(r"D:\model_research\yeast_GEM_multi_omics_analysis\N_lim\output\p distri1.tif", dpi=600)
    # plt.savefig('./output/p distri1.tif', dpi=600)

    plt.show()
