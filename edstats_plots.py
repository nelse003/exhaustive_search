import pandas as pd
import matplotlib
matplotlib.use('agg')
import seaborn as sns

es_minima_csv = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus/es_minima.csv"
edstats_csv = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus/edstats.csv"

edstats_df = pd.read_csv(edstats_csv)
es_minima_df = pd.read_csv(es_minima_csv, names=['Dataset','ES_occ','es_bfac','min_fofc'])

summary_df = pd.merge(edstats_df,es_minima_df, on='Dataset')

print(edstats_df)
print(es_minima_df)

print(edstats_df['Dataset'])
print(es_minima_df['Dataset'])

print(pd.merge(edstats_df,es_minima_df, on='Dataset'))
print(list(pd.merge(edstats_df,es_minima_df, on='Dataset')))

pairplot = sns.pairplot(summary_df)
sns_plot.savefig("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus/pairplot.png",
                 dpi=600)