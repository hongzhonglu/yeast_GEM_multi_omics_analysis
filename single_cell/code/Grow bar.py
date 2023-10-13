##################################################
# get "Grow violin.jpg"
# first run "get_pfba_flux.py" to compute flux data "pfba_flux.xlsx"

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os


data = pd.read_excel('../output/pfba_flux.xlsx')
stress = data.loc[0:77, 'r_2111']
unstress = data.loc[78:, 'r_2111']

# Dataset:
a = pd.DataFrame({'group': np.repeat('Stress', 80), 'value': data.loc[0:79, 'r_2111']})
b = pd.DataFrame({'group': np.repeat('Unstress', 83), 'value': data.loc[80:, 'r_2111']})
df = b.append(a)
# draw
font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 23,
         }
font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 24,
         }
font3 = {'family': 'Arial',
         'weight': 'normal',
         'size': 24,
         }
# plot violin chart
plt.subplots(figsize=(8, 6))
sns.violinplot(x='group',
               y='value',
               data=df,
               )
# add title
#plt.ylim(0.08, 0.0852)
plt.xticks(fontproperties='Arial', size=24)
plt.yticks(fontproperties='Arial', size=20)
plt.xlabel('')
plt.ylabel('Growth rate(h-1)',
           fontdict=font3)
plt.tight_layout()
plt.savefig('../output/Grow violin.jpg',
            dpi=800)
plt.show()
