import os
import pandas as pd
import matplotlib
import sqlite3
matplotlib.use('agg')
import seaborn as sns
from matplotlib import pyplot as plt

def labelled_pairplot(df, hue_column=None):

    if hue_column is not None:
        g = sns.PairGrid(df, hue = hue_column)
    else:
        g = sns.PairGrid(df)

    g = g.map_diag(plt.hist)
    g = g.map_offdiag(plt.scatter)

    xlabels, ylabels = [], []

    for ax in g.axes[-1, :]:
        xlabel = ax.xaxis.get_label_text()
        xlabels.append(xlabel)
    for ax in g.axes[:, 0]:
        ylabel = ax.yaxis.get_label_text()
        ylabels.append(ylabel)

    for i in range(len(xlabels)):
        for j in range(len(ylabels)):
            g.axes[j, i].xaxis.set_label_text(xlabels[i])
            g.axes[j, i].yaxis.set_label_text(ylabels[j])

    return g


es_minima_csv = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus/es_minima.csv"
edstats_csv = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus/edstats.csv"
out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus"
database_path = "/dls/labxchem/data/2016/lb13385-64/processing/database/soakDBDataFile.sqlite"

# out_dir =  "/home/nelse003/Desktop/DCP2B_exhaus/"
#
# es_minima_csv = "/home/nelse003/Desktop/DCP2B_exhaus/es_minima.csv"
# edstats_csv = "/home/nelse003/Desktop/DCP2B_exhaus/edstats.csv"
#database_path = "/home/nelse003/Desktop/DCP2B_exhaus/soakDBDataFile.sqlite"

edstats_df = pd.read_csv(edstats_csv)
es_minima_df = pd.read_csv(es_minima_csv, names=['Dataset','ES_occ','es_bfac','min_fofc'])

summary_df = pd.merge(edstats_df,es_minima_df, on='Dataset')

refinement_outcomes= "'4 - CompChem ready', '5 - Deposition ready','6 - Deposited'"

print(database_path)
conn = sqlite3.connect(database_path)
main_table_df = pd.read_sql_query("select * from mainTable",conn)
cur = conn.cursor()

cur.execute("SELECT CrystalName, CompoundCode, RefinementResolution "
            "FROM mainTable WHERE RefinementOutcome in ({})"
            " AND  (RefinementPDB_latest AND RefinementMTZ_latest) IS NOT NULL".format(refinement_outcomes))

refinement_xtals = cur.fetchall()

#Close connection to the database
cur.close()

summary_df = summary_df.rename(columns={'Dataset':'CrystalName'})
# summary_df = pd.merge(summary_df, main_table_df, on='CrystalName')


compounds = {}
for xtal_name, compound_code, resolution in refinement_xtals:
    #xtal_name = xtal_name.encode('ascii')
    #compound_code = compound_code.encode('ascii')
    compounds[xtal_name] = compound_code

comp_df = pd.DataFrame(list(compounds.items()), columns=['CrystalName','compound_code'])
summary_df = pd.merge(summary_df, comp_df, on='CrystalName')

FMOPL000435a_df = summary_df[summary_df['compound_code'] == "FMOPL000435a"]

if not os.path.exists(os.path.join(out_dir,"pairplot.png")):
    pairplot = labelled_pairplot(summary_df, hue_column="compound_code")
    fig = pairplot.fig
    fig.savefig(os.path.join(out_dir,"pairplot.png"),
                     dpi=300)

if not os.path.exists(os.path.join(out_dir,"FMOPL000435a_pairplot.png")):
    FMOPL000435a_pairplot = labelled_pairplot(FMOPL000435a_df)
    fig = FMOPL000435a_pairplot.fig
    fig.savefig(os.path.join(out_dir,"FMOPL000435a_pairplot.png"), dpi=300)

summary_df.to_csv(os.path.join(out_dir, "DCP2B_edstats_summary.csv"))
FMOPL000435a_df.to_csv(os.path.join(out_dir, "FMOPL000435a_edstats_summary.csv"))

duplicate_compound_df = pd.concat(g for _, g in summary_df.groupby("compound_code") if len(g) > 1)

summary_duplicate_df_list = []
for duplicate_compound in duplicate_compound_df['compound_code'].unique():
    duplicate_df =  summary_df[summary_df['compound_code'] == duplicate_compound]

    summary = {"compound" : [duplicate_compound],
               "number refined hits" : [len(duplicate_df.index)],
               "RSCC min" : [duplicate_df['RSCC'].min()],
               "RSCC max" : [duplicate_df['RSCC'].max()],
               "Occ refined min" : [duplicate_df['Occupancy'].min()],
               "Occ refined max" : [duplicate_df['Occupancy'].max()],
               "Occ ES min": [duplicate_df['ES_occ'].min()],
               "Occ ES max": [duplicate_df['ES_occ'].max()]
               }
    summary_duplicate_df_list.append(pd.DataFrame(data=summary))

    # print("{} length: {}".format(duplicate_compound, len(duplicate_df.index)))
    # print("{} RSCC range: {} {}").format(duplicate_compound, duplicate_df['RSCC'].min(),
    #         duplicate_df['RSCC'].max())
    # print("{} refinement_occ range: {} {}").format(duplicate_compound,
    #         duplicate_df['Occupancy'].min(), duplicate_df['Occupancy'].max())
    # print("{} es_occ range: {} {}").format(duplicate_compound,
    #         duplicate_df['ES_occ'].min(), duplicate_df['ES_occ'].max())
    #
    # print("_______________________________________________________________________")

pd.concat(summary_duplicate_df_list).to_csv(os.path.join(out_dir, "DCP2B_edstats_duplicates.csv"))