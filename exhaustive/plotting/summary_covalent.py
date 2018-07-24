import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

cvr = pd.read_csv('/home/nelse003/Documents/covalent_ratio_subjective.csv',
                  names=['ratio','cateogry'], engine='python')

ct = pd.crosstab(index=cvr["ratio"],columns=cvr["cateogry"])

colors=['xkcd:grey','xkcd:green', 'xkcd:orange']

ct.plot.bar(stacked=True, color=colors, legend=False)
plt.xlabel("Percentage of labelled protein")
plt.ylabel("Number of Crystals")
plt.legend(loc='best')
plt.title("Strength of 2mFo-DFc evidence for covalent ligand")

# stacked = ct.stack().reset_index().rename(columns={0:'value'})
# sns_plot = sns.barplot(x=stacked.cateogry, y=stacked.value, hue=stacked.mark)

plt.savefig('/home/nelse003/Documents/example.png')
