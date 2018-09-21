import pandas as pd

es_minima_csv = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus/es_minima.csv"
edstats_csv = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus/edstats.csv"

edstats_df = pd.read_csv(edstats_csv)
es_minima_df = pd.read_csv(es_minima_csv, names=['Dataset','ES_occ','es_bfac''min_fofc'])

print(edstats_df)
print(es_minima_df)

print(pd.concat([edstats_df,es_minima_df], axis=0))

