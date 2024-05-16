import numpy as np
import matplotlib.pyplot as plt

# width of the bars
barWidth = 0.3

# Choose the height of the blue bars
# bars1 = [88.09, 66]
#
# # Choose the height of the cyan bars
# bars2 = [85.75, 52]
bars1 = [88.09, 66]

# Choose the height of the cyan bars
bars2 = [85.75, 52]

# The x position of bars
r1 = np.arange(len(bars1))
r2 = [x + 1.2*barWidth for x in r1]

# Create blue bars
plt.bar(r2, bars1, width=barWidth, color='#6f9bc6', edgecolor='black', capsize=7, label='Yeast9')

# Create cyan bars
plt.bar(r1, bars2, width=barWidth, color='#f3993a', edgecolor='black', capsize=7, label='Yeast8.3.0')

font1 = {'family': 'Arial',
         'weight': 'normal',
         'size': 20,
         }
font2 = {'family': 'Arial',
         'weight': 'normal',
         'size': 24,
         }

# general layout
plt.xticks([r + 1.2*barWidth/2 for r in range(len(bars1))],
           ['Gene essentiality', 'Memote total score'],
           )
plt.ylabel('Percentage (%)', fontdict=font2)
plt.yticks(np.arange(0, 100, 10), fontproperties='Arial', size=20)
plt.xticks(fontproperties='Arial', size=22)

font = {'family': 'Arial', 'weight': 'normal', 'size': 20}
legend = plt.legend(prop=font)
# Show graphic
plt.tight_layout()
plt.savefig('./memote.tif', dpi=600)
plt.show()
