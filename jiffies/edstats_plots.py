import os
import sqlite3

import matplotlib
import pandas as pd

matplotlib.use('agg')
import seaborn as sns
from matplotlib import pyplot as plt
from utils import u_iso_to_b_fac, get_minimum_fofc
from exhaustive import master_phil
from plot import occupancy_histogram_with_exhaustive_search


def labelled_pairplot(df, hue_column=None):
    """ Seaborn pairplot with labelled axis for each subplot
    
    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe to be plotted 
    
    hue_column: str
        Passed to sns.pairgrid to set the colour of points based on 
        categorisation by values in a column of the dataframe
        
    Returns
    -------
    g: seaborn.pairGrid
        pair grid object that can then be plotted

    """

    if hue_column is not None:
        g = sns.PairGrid(df, hue=hue_column)
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


def plot_NUDT22():
    NUDT22_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/NUDT22_repeat_soaks"
    prefix = "NUDT22A-x"

    compound_dirs = [os.path.join(NUDT22_dir, compound_dir) for compound_dir in os.listdir(NUDT22_dir)
                     if os.path.isdir(os.path.join(NUDT22_dir, compound_dir))]

    compound_dirs.remove(os.path.join(loop_dir, 'N14004a'))

    combined_occ_df_list = []
    for compound_dir in compound_dirs:
        pass


def plot_DCP2B():
    """ Plotting code for edtstats pairplots on DCP2B 
    
    Notes
    -------------
    Place holder to be better functioanlised when 
    more generalised case is apparent
    
    """

    es_minima_csv = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus/es_minima.csv"
    edstats_csv = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus/edstats.csv"
    out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus"
    database_path = "/dls/labxchem/data/2016/lb13385-64/processing/database/soakDBDataFile.sqlite"

    edstats_df = pd.read_csv(edstats_csv)
    es_minima_df = pd.read_csv(es_minima_csv, names=['Dataset', 'ES_occ', 'es_bfac', 'min_fofc'])

    summary_df = pd.merge(edstats_df, es_minima_df, on='Dataset')

    refinement_outcomes = "'4 - CompChem ready', '5 - Deposition ready','6 - Deposited'"

    print(database_path)
    conn = sqlite3.connect(database_path)
    main_table_df = pd.read_sql_query("select * from mainTable", conn)
    cur = conn.cursor()

    cur.execute("SELECT CrystalName, CompoundCode, RefinementResolution "
                "FROM mainTable WHERE RefinementOutcome in ({})"
                " AND  (RefinementPDB_latest AND RefinementMTZ_latest) IS NOT NULL".format(refinement_outcomes))

    refinement_xtals = cur.fetchall()

    # Close connection to the database
    cur.close()

    summary_df = summary_df.rename(columns={'Dataset': 'CrystalName'})
    # summary_df = pd.merge(summary_df, main_table_df, on='CrystalName')

    compounds = {}
    es_occ_b = []
    for xtal_name, compound_code, resolution in refinement_xtals:
        # xtal_name = xtal_name.encode('ascii')
        # compound_code = compound_code.encode('ascii')
        compounds[xtal_name] = compound_code

        if compound_code == "FMOPL000435a":
            csv = os.path.join(out_dir, xtal_name, "exhaustive_search.csv")
            occ, u_iso, _ = get_minimum_fofc(csv)
            es_occ_b.append([xtal_name,
                             occ,
                             u_iso_to_b_fac(u_iso)])

    comp_df = pd.DataFrame(list(compounds.items()), columns=['CrystalName', 'compound_code'])
    summary_df = pd.merge(summary_df, comp_df, on='CrystalName')

    FMOPL000435a_df = summary_df[summary_df['compound_code'] == "FMOPL000435a"]

    if not os.path.exists(os.path.join(out_dir, "pairplot.png")):
        pairplot = labelled_pairplot(summary_df, hue_column="compound_code")
        fig = pairplot.fig
        fig.savefig(os.path.join(out_dir, "pairplot.png"),
                    dpi=300)

    if not os.path.exists(os.path.join(out_dir, "FMOPL000435a_pairplot.png")):
        FMOPL000435a_pairplot = labelled_pairplot(FMOPL000435a_df)
        fig = FMOPL000435a_pairplot.fig
        fig.savefig(os.path.join(out_dir, "FMOPL000435a_pairplot.png"), dpi=300)

    FMOPL000435a_df = FMOPL000435a_df.rename(index=str, columns={"ES_occ": "es_occupancy",
                                                                 "Occupancy": "occupancy"})

    params = master_phil.extract()
    params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus"

    occupancy_histogram_with_exhaustive_search(FMOPL000435a_df,
                                               protein_name="DCP2B",
                                               compound="FMOPL000435a",
                                               params=params)
    # occupancy_b_factor_scatter_plot(FMOPL000435a_df,
    #                                 protein_name=protein_name,
    #                                 compound=compound,
    #                                 params=params)

    summary_df.to_csv(os.path.join(out_dir, "DCP2B_edstats_summary.csv"))
    FMOPL000435a_df.to_csv(os.path.join(out_dir, "FMOPL000435a_edstats_summary.csv"))
    duplicate_compound_df = pd.concat(g for _, g in summary_df.groupby("compound_code") if len(g) > 1)

    summary_duplicate_df_list = []
    for duplicate_compound in duplicate_compound_df['compound_code'].unique():
        duplicate_df = summary_df[summary_df['compound_code'] == duplicate_compound]

        summary = {"compound": [duplicate_compound],
                   "number refined hits": [len(duplicate_df.index)],
                   "RSCC min": [duplicate_df['RSCC'].min()],
                   "RSCC max": [duplicate_df['RSCC'].max()],
                   "Occ refined min": [duplicate_df['Occupancy'].min()],
                   "Occ refined max": [duplicate_df['Occupancy'].max()],
                   "Occ ES min": [duplicate_df['ES_occ'].min()],
                   "Occ ES max": [duplicate_df['ES_occ'].max()]
                   }
        summary_duplicate_df_list.append(pd.DataFrame(data=summary))

    pd.concat(summary_duplicate_df_list).to_csv(os.path.join(out_dir, "DCP2B_edstats_duplicates.csv"))


plot_DCP2B()
