import argparse
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd

def u_iso_to_b_fac(u_iso):
    """
    Convert isotropic B-factor to

    Parameters
    ----------
    u_iso: float
        isotropic displacement (u_iso)

    Returns
    -------
    b_iso: float
        isotropic B-factor
    """
    b_iso = (8 * np.pi ** 2) * u_iso
    return b_iso


def scatter_plot(
    csv_name, three_dim_plot=True, occ_plot=False, b_plot=False, title_text=None
):

    """
    Scatter plots of occupancy, U_iso and mean |fo_fc| from csv output

    Parameters
    ----------
    csv_name: str
        path to csv file

    three_dim_plot: Bool
        flag to determine whether a three dimensional plot is used

    title_text:str
        plot title text

    Returns
    -------
    None

    Notes
    ------
    Saves png to same name as csv
    """
    # TODO Split into multiple functions with base class

    # Load data from CSV
    if csv_name.endswith(".csv"):
        data = np.genfromtxt(csv_name, delimiter=",", skip_header=0)
    else:
        data = np.genfromtxt("{}.csv".format(csv_name), delimiter=",", skip_header=0)

    if len(data[0]) == 3:
        occ = data[:, 0]
        u_iso = data[:, 1]
        fo_fc = data[:, 2]

    if len(data[0]) == 4:
        occ = data[:, 0]
        u_iso = data[:, 2]
        fo_fc = data[:, 3]

    b_iso = u_iso_to_b_fac(u_iso)

    fig = plt.figure()

    if three_dim_plot:
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(occ, b_iso, fo_fc)
        plt.xlabel("Occupancy")
        plt.ylabel("B Factor")
        ax.set_zlabel("Mean(|mFo-Fc|)")
        type = "3d"

    elif occ_plot:
        ax = fig.add_subplot(111)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.scatter(occ, fo_fc)

        plt.xlabel("Occupancy")
        plt.ylabel("Mean(|mFo-DFc|)")
        type = "occ"

    elif b_plot:
        ax = fig.add_subplot(111)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.scatter(b_iso, fo_fc)

        cbar.ax.set_ylabel('# of contacts', rotation=270)
        type = "b"

    if title_text is not None:
        plt.title(title_text)

    plt.savefig("{}_{}".format(csv_name.rstrip(".csv"), type), dpi=300)
    plt.close()

def color_occ_scatter_plot(csv_name):

    df = pd.read_csv(csv_name,
                     header=None,
                     names=["bound_occupancy",
                            "ground_occupancy",
                            "u_iso",
                            "fo_fc"])

    plt.scatter(x=u_iso_to_b_fac(df.u_iso),
                y=df.fo_fc,
                c=df.bound_occupancy,
                cmap="Blues")
    cbar =plt.colorbar()
    cbar.ax.set_ylabel('Occupancy', rotation=270)
    cbar.ax.get_yaxis().labelpad = 15
    plt.xlabel("B factor")
    plt.ylabel("Mean(|mFo-DFc|)")
    plt.savefig("occ_colour_DCP2B-x0146", dpi=300)

if __name__ == "__main__":
    """Create plots from exhaustive csv file
    
    This file was not oworking in ccp4-python due to matplotlibe version issues.
    Added functions required to run in exhaustive conda environemnt instead, 
    and removed depedency on other scripts, by copying from utils.utils
    """

    # Add ability to parse command line arguments
    parser = argparse.ArgumentParser(
        description="Create plots from exhaustive csv file"
    )
    # Add csv path argument
    parser.add_argument("--csv", action="store", dest="csv")

    # resolve the parser
    args = parser.parse_args()

    color_occ_scatter_plot(csv_name=args.csv)
    exit()

    scatter_plot(csv_name=args.csv)
    scatter_plot(csv_name=args.csv, occ_plot=True, three_dim_plot=False)
    scatter_plot(csv_name=args.csv, b_plot=True, three_dim_plot=False)
