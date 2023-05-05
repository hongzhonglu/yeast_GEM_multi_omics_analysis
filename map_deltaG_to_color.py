import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def num2color(values, cmap):
    """Map values to colors"""
    norm = mpl.colors.Normalize(vmin=np.min(values), vmax=np.max(values))
    cmap = mpl.cm.get_cmap(cmap)
    return [cmap(norm(val)) for val in values]

def to_RGB(color):
    return color * 255

big = pd.read_table('morethan0.txt')
small = pd.read_table('lessthan0.txt')
small.loc[5, 'g'] = 0.0

big.loc[:, 'g'] = -big.loc[:, 'g']
bigcolor = num2color(big.loc[:, 'g'], "Blues")
bigcolor = pd.DataFrame(bigcolor, index=big.loc[:, 'id'], columns=['R', 'G', 'B', 'T'])
bigcolor = bigcolor.apply(to_RGB)
bigcolor.to_excel('bigcolorRGB.xlsx')

# small.loc[:, 'g'] = -small.loc[:, 'g']
small.loc[:, 'g'] = -small.loc[:, 'g']
smallcolor = num2color(small.loc[:, 'g'], "Reds")
smallcolor = pd.DataFrame(smallcolor, index=small.loc[:, 'id'], columns=['R', 'G', 'B', 'T'])
smallcolor = smallcolor.apply(to_RGB)
smallcolor.to_excel('smallolorRGB.xlsx')



