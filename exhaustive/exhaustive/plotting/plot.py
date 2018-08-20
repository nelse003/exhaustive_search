from __future__ import division, print_function
import os
import matplotlib
import numpy as np
import logging
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from exhaustive.exhaustive.utils.utils import (get_fofc_from_csv,
get_minimum_fofc, round_step, b_to_u_iso, u_iso_to_b_fac)
from mpl_toolkits.mplot3d import Axes3D
#import seaborn as sns
##############################################################

logging = logging.getLogger(__name__)


def scatter_plot(csv_name, three_dim_plot=True, title_text=None):

    """ Scatter plots of occupancy, U_iso and mean |fo_fc| from csv output """

    # Load data from CSV
    if csv_name.endswith(".csv"):
        data = np.genfromtxt(csv_name, delimiter=',', skip_header=0)
    else:
        data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',',
         skip_header=0)

    if len(data[0]) == 3:
        logging.info("Using 3 column data to plot")
        occ = data[:, 0]
        u_iso = data[:, 1]
        fo_fc = data[:, 2]

    if len(data[0]) == 4:
        logging.info("Using 4 column data to plot")
        occ = data[:, 0]
        u_iso = data[:, 2]
        fo_fc = data[:, 3]

    b_iso = u_iso_to_b_fac(u_iso)

    fig = plt.figure()

    if three_dim_plot:
        logging.info("Plotting 3d plot")
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(occ, b_iso, fo_fc)
        plt.xlabel("Occupancy")
        plt.ylabel("B_iso")
        ax.set_zlabel("Fo-Fc")
    else:
        ax = fig.add_subplot(111)
        ax.scatter(occ, fo_fc)

        plt.xlabel("Occupancy")
        plt.ylabel("Fo-Fc")

    if title_text is not None:
        plt.title(title_text)

    # TODO Plot name parameter #54, requires templating
    plt.savefig(csv_name.rstrip(".csv"))
    plt.close()


def bounded_2d_scatter(atom_name, lower_bound, upper_bound):

    # Load data from CSV
    data = np.genfromtxt('{}.csv'.format(atom_name), delimiter=',',
                         skip_header=0)

    reduced_data = data[np.where(np.logical_and(data[:, 2] >= lower_bound,
    data[:, 2] <= upper_bound))]

    x = reduced_data[:, 0]
    y = reduced_data[:, 1]

    y = (8 * np.pi ** 2) * y ** 2

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlabel("Occupancy")
    plt.ylabel("B_iso")
    plt.title("Fofc values between {} and {}".format(lower_bound, upper_bound))

    ax.scatter(x, y)
    plt.savefig("{}_reduced_".format(atom_name))
    plt.close()


def colourbar_2d_scatter(atom_name):
    # Load data from CSV
    data = np.genfromtxt('{}.csv'.format(atom_name),
                         delimiter=',', skip_header=0)

    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    y = (8 * np.pi ** 2) * y ** 2

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlabel("Occupancy")
    plt.ylabel("B_iso")

    sc = ax.scatter(x, y, c=z, marker=".")

    plt.colorbar(sc)
    plt.savefig("{}_colorbar".format(atom_name))
    plt.close()


def connectpoints(x, y, x_1, y_1, p1, linestyle='k-'):

    """ Draw lines between two lists of points"""

    x1, x2 = x[p1], x_1[p1]
    y1, y2 = y[p1], y_1[p1]
    plt.plot([x1, x2], [y1, y2], linestyle)


def connectpoint(x, y, x_1, y_1, linestyle='k-'):

    """ Draw line between two points """
    plt.plot([x, x_1], [y, y_1], linestyle)


def connectpoints_3d(x, y, z, x_1, y_1, z_1, p1, linestyle='k-'):

    """ Draw lines between two sets od points in 3d"""

    x1, x2 = x[p1], x_1[p1]
    y1, y2 = y[p1], y_1[p1]
    z1, z2 = z[p1], z_1[p1]
    plt.plot([x1, x2], [y1, y2], [z1, z2], linestyle)


# TODO sort out params
def plot_3d_fofc_occ(start_occ,
                     end_occ,
                     step,
                     dataset_prefix,
                     set_b,
                     out_dir,
                     params):

    """ Plot the difference in occupancy & mean(|fo-fc|)
    at the simulated occupancy and the minima. """

    min_fofcs = []
    min_occs = []
    min_b_facs = []
    fofcs = []
    occs = []
    b_facs = []

    for lig_occupancy in np.arange(start_occ, end_occ + (step / 5), step):

        # TODO Replace CSV naming #59

        csv_name = params.exhaustive.output.csv_prefix \
                   + "_occ_{}_b_{}.csv".format(
                   str(lig_occupancy).replace(".", "_"),
                   str(set_b).replace(".", "_"))

        csv_path = os.path.join(out_dir, csv_name)

        min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_path)
        fofc = get_fofc_from_csv(csv_name, lig_occupancy,
                                 round_step(b_to_u_iso(set_b)),
                                 step)
        fofcs.append(fofc)
        occs.append(lig_occupancy)
        b_facs.append(set_b)
        min_b_facs.append(u_iso_to_b_fac(min_u_iso))
        min_fofcs.append(fo_fc_at_min)
        min_occs.append(min_occ)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    min_plot, = ax.plot(min_occs, min_b_facs, min_fofcs, 'k+')
    occ_plot, = ax.plot(occs, b_facs, fofcs, 'ro')

    for i in np.arange(0, len(occs)):
        connectpoints_3d(occs, b_facs, fofcs, min_occs, min_b_facs,
                         min_fofcs, i)

    ax.legend((min_plot, occ_plot),
              ('Minima of mean |Fo-Fc|',
               'Mean |Fo-Fc| at simulated occupancy'),
              prop={"size": 8},
              numpoints=1,
              bbox_to_anchor=(1, 1),
              bbox_transform=plt.gcf().transFigure)

    ax.set_xlabel("Occupancy")
    ax.set_ylabel("B factor")
    ax.set_zlabel("Mean |Fo-Fc|")

    plt.title("{}: Delta mean|Fo-Fc| and Delta Occupancy".format(
        dataset_prefix), fontsize=10)

    # TODO Replace with vlaidate.output.plot_name, requires templating #54
    plt.savefig("{}-3d-delta_fofc_occ.png".format(dataset_prefix))


def occupancy_histogram_with_exhaustive_search(occ_df, protein_name, compound,
                                               params):

    es_occs = occ_df['es_occupancy'].dropna()
    refined_occs = occ_df['occupancy'].dropna()

    if len(es_occs) == 0 | len(refined_occs) == 0:
        print("No results to plot: {} {}".format(protein_name, compound))
        return

    print(refined_occs)
    print(es_occs)

    occ_bins = np.linspace(0, 1, 21, endpoint=True)

    fig, ax = plt.subplots()
    try:
        ax.hist(es_occs, bins=occ_bins, width=0.04, color='r',  alpha=0.5,
                label='Exhaustive search occupancy: {}'.format(len(es_occs)))
        ax.hist(refined_occs, bins=occ_bins, width=0.04, color='b', alpha=0.5,
                label='Refined occupancy: {}'.format(len(refined_occs)))
    except ValueError:
        print("Can't make histogram for: {} {}".format(protein_name, compound))
        return

    plt.xlim(0, 1.05)
    ax.legend(loc='best', fontsize='small')

    plt.title("Occupancy Histogram: {} : {}".format(protein_name, compound))
    plt.xlabel("Occupancy")
    plt.ylabel("Frequency")

    file_path = os.path.join(params.output.out_dir,
                             "occ_hist_with_exhaustive_{}_{}.png".format(
                                 protein_name, compound))

    plt.savefig(file_path, dpi=300)
    plt.close()


def occupancy_b_factor_scatter_plot(occ_df, protein_name, compound, params):

    es_occs = occ_df['es_occupancy']
    refined_occs = occ_df['occupancy']
    es_b_fac = occ_df['es_b_fac']
    refine_mean_b_fac = occ_df['mean_b_fac']
    refine_std_b_fac = occ_df['std_b_fac']

    es_occs_joined = occ_df.dropna()['es_occupancy']
    refined_occs_joined = occ_df.dropna()['occupancy']
    es_b_fac_joined = occ_df.dropna()['es_b_fac']
    refine_mean_b_fac_joined = occ_df.dropna()['mean_b_fac']

    fig, ax = plt.subplots()

    ax.scatter(es_occs, es_b_fac,
               label="Exhaustive search: {}".format(
                   len(occ_df['es_occupancy'].dropna())))

    ax.errorbar(refined_occs, refine_mean_b_fac,
                fmt='rs',
                yerr=refine_std_b_fac,
                label="Refinement (errorbar = standard deviation of B factor "
                     "across ligand): {}".format(
                    len(occ_df['occupancy'].dropna())),
                linestyle="None")

    for i in np.arange(0, len(es_occs_joined)):
        connectpoints(es_occs_joined, es_b_fac_joined, refined_occs_joined,
                      refine_mean_b_fac_joined, i, linestyle='k:')

    # TODO Add post Exhaustive search refinement

    ax.legend(loc='best', fontsize='small')
    plt.xlabel("Occupancy")
    plt.ylabel("B Factor")
    plt.title("Exhaustive search minima compared to refined "
              "{} : {}".format(protein_name, compound))

    file_path = os.path.join(params.output.out_dir,
                             "occ_scatter_with_exhaustive_"
                             "{}_{}".format(protein_name,  compound))

    plt.savefig(file_path, dpi=300)
    plt.close()

    print(occ_df)

    # if len(occ_df) != len(occ_df.dropna()):
    #     print("!!!!!!")
    #     print(occ_df)
    #     print(occ_df.dropna())
    #     exit()

def plot_fofc_occ(start_occ, end_occ, step, dataset_prefix, set_b):

    """Plot the difference in occupancy/fofc at the simulated occupancy
    and minima."""

    min_fofcs = []
    min_occs = []
    fofcs = []
    occs = []

    for lig_occupancy in np.arange(start_occ, end_occ + (step / 5), step):

        csv_name = "occ_{}_b_{}_u_iso".format(
            str(lig_occupancy).replace(".", "_"),
            str(set_b).replace(".", "_"))

        min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_name)

        fofc = get_fofc_from_csv(csv_name,
                                 lig_occupancy,
                                 round_step(b_to_u_iso(set_b)),
                                 step)
        fofcs.append(fofc)
        occs.append(lig_occupancy)
        min_fofcs.append(fo_fc_at_min)
        min_occs.append(min_occ)

    fig, ax = plt.subplots()
    min_plot, = ax.plot(min_occs, min_fofcs, 'k+')
    occ_plot, = ax.plot(occs, fofcs, 'ro')

    for i in np.arange(0, len(occs)):
        connectpoints(occs, fofcs, min_occs, min_fofcs, i)

    ax.legend((min_plot, occ_plot),
              ('Minima of mean |Fo-Fc|',
               'Mean |Fo-Fc| at simulated occupancy'),
              prop={"size": 8},
              numpoints=1,
              bbox_to_anchor=(1, 1),
              bbox_transform=plt.gcf().transFigure)

    ax.set_xlabel("Occupancy")
    ax.set_ylabel("Mean |Fo-Fc|")

    plt.title("{}: Delta mean|Fo-Fc| "
              "and Delta Occupancy".format(dataset_prefix), fontsize=10)
    plt.savefig("{}-delta_fofc_occ.png".format(dataset_prefix))

def plot_edstats_metric(edstats_df, compound_folder, compound, protein_name,
                        metric_name):

    """Plot a single giant.score_model metric"""

    metric_filename = metric_name.replace(" ", "_")
    metric_filename = metric_filename.replace("/", "_")

    filename = os.path.join(
        compound_folder,
        "{}_{}_{}.png".format(protein_name, compound, metric_filename))

    datasets = edstats_df['Dataset'].tolist()

    ax=edstats_df.plot(x="Dataset", y=metric_name, marker='o', linestyle='None')
    ax.legend().set_visible(False)
    plt.title("{} {}".format(protein_name, compound))
    plt.ylabel(metric_name)
    plt.xticks(np.arange(len(datasets)),datasets,rotation='vertical')
    plt.xlim(-1,len(datasets)+1)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_edstats_across_soaks(edstats_df, compound_folder, compound,
                              protein_name, title_suffix=None):

    metrics = ["Occupancy", "RSZO/OCC", "Model RMSD", "RSZD", "RSR", "RSCC",
                   "Surroundings B-factor Ratio"]
    metrics_2 = ["Occupancy-2", "RSZO/OCC-2", "RSZD-2", "RSR-2", "RSCC-2"]

    for metric in metrics + metrics_2:
        plot_edstats_metric(edstats_df, compound_folder, compound, protein_name,
                            metric)

    #pairplot



    cols = [c for c in edstats_df.columns if c.lower()[-2:] == '-2']
    cols_1 = [col.rstrip('-2') for col in cols]
    metric_df = edstats_df.drop(cols, axis=1)
    metric_2_df = edstats_df.drop(cols_1, axis=1)

    metric_df['source_file'] = 1
    metric_2_df['source_file'] = 2

    for i in np.arange(len(cols)):
        metric_2_df = metric_2_df.rename(columns={cols[i]:cols_1[i]})

    edstats_df = pd.concat([metric_df, metric_2_df], ignore_index=True)

    sns.pairplot(edstats_df, hue="source_file")
    plt.tight_layout()
    plt.savefig(os.path.join(compound_folder, "pairplot.png"))

    pp = sns.pairplot(edstats_df, vars=metrics, hue="source_file")
    pp.fig.suptitle("{} {} :{}".format(protein_name,compound,title_suffix))
    plt.tight_layout()
    plt.savefig(os.path.join(compound_folder, "reduced_pairplot.png"))




