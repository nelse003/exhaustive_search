import matplotlib
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt

cvr = pd.read_csv('/home/nelse003/Documents/covalent_ratio_subjective.csv',
                  names=['ratio', 'cateogry'], engine='python')

ct = pd.crosstab(index=cvr["ratio"], columns=cvr["cateogry"])

colors = ['xkcd:grey', 'xkcd:green', 'xkcd:orange']

ct.plot.bar(stacked=True, color=colors, legend=False)
plt.xlabel("Percentage of labelled protein")
plt.ylabel("Number of Crystals")
plt.legend(loc='best')
plt.title("Strength of 2mFo-DFc evidence for covalent ligand")
plt.tight_layout()

plt.savefig('/home/nelse003/Documents/example.png', dpi=600, transparent=True)
