##################################################
# get "Grow violin.jpg"
# first run "get_pfba_flux.py" to compute flux data "pfba_flux.xlsx"

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

os.chdir('..')

data = pd.read_excel('../singlecell/output/pfba_flux.xlsx')
stress = data.loc[0:79, 'r_2111']
unstress = data.loc[80:, 'r_2111']

# Dataset:
a = pd.DataFrame({'group': np.repeat('stress', 80), 'value': data.loc[0:79, 'r_2111']})
b = pd.DataFrame({'group': np.repeat('unstress', 83), 'value': data.loc[80:, 'r_2111']})
df = a.append(b)
# draw
font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 23,
         }
font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 20,
         }
font3 = {'family': 'Arial',
         'weight': 'normal',
         'size': 16,
         }
# plot violin chart
sns.violinplot(x='group',
               y='value',
               data=df,
               )
# add title
plt.xticks(fontproperties='Arial', size=16)
plt.yticks(fontproperties='Arial', size=16)
plt.xlabel('')
plt.ylabel('Growth rate(h-1)',fontdict=font3)
plt.tight_layout()
plt.savefig('../singlecell/output/Grow violin.jpg',
            dpi=600)
plt.show()
