import os

import numpy as np
import pandas as pd


def b_to_u_iso(b_fac):
    """
    Convert isotropic displacement (u_iso) to isotropic B-factor

    Parameters
    ----------
    b_fac: float
        isotropic B-factor

    Returns
    -------
    u_iso: float

    """
    u_iso = b_fac / (8 * np.pi ** 2)
    return u_iso


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


def get_minimum_fofc(csv_name, b_fac=None):
    """
    Get minima in fofc, and return minima and where it occurs

    B factor can be supplied to look at the minima across a
    single b factor value

    Parameters
    ----------

    Returns
    --------
    """
    # TODO Remove this dual .csv by cleaning up csv name: issue 59
    if csv_name.endswith(".csv"):
        data = np.genfromtxt(csv_name, delimiter=",", skip_header=0)
    else:
        data = np.genfromtxt("{}.csv".format(csv_name), delimiter=",", skip_header=0)

    # If four column data from multiple ligand
    if len(data[0]) == 4:
        occ = data[:, 0]
        u_iso = data[:, 2]
        fo_fc = data[:, 3]
    elif len(data[0]) == 3:
        occ = data[:, 0]
        u_iso = data[:, 1]
        fo_fc = data[:, 2]
    else:
        raise ValueError("|Fo-Fc| data is not in correct format")

    if b_fac is not None:
        set_u_iso = b_to_u_iso(b_fac)
        step = np.unique(occ)[1] - np.unique(occ)[0]
        check_u_iso = round_step(set_u_iso, base=step)
        data_array = np.stack((occ, u_iso, fo_fc))

        occ = data_array[0][
            (data_array[1] >= check_u_iso) & (data_array[1] <= check_u_iso)
        ]

        u_iso = data_array[1][
            (data_array[1] >= check_u_iso) & (data_array[1] <= check_u_iso)
        ]

        fo_fc = data_array[2][
            (data_array[1] >= check_u_iso) & (data_array[1] <= check_u_iso)
        ]

    min_index = np.argmin(fo_fc)

    return occ[min_index], u_iso[min_index], fo_fc[min_index]


def round_step(x, prec=2, base=0.05):
    """ Return a number rounded to the nearest base."""
    return round(base * round(float(x) / base), prec)


def expand_array(array):
    "Expand a 3d numpy array"

    x = array[:, 0]
    y = array[:, 1]
    z = array[:, 2]

    return x, y, z


def sample_spherical(npoints, ndim=3):
    """Sample a ndimensional sphere using gaussians

    Is used to generate points within a sphere from which
    a convex hull around atoms can be generated.

    Parameters
    ----------
    npoints:

    ndim: int

    Returns
    -------
    """

    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec


def wait_for_file_existence(file_path, wait_time):
    """ Wait for a file to exist, stop after wait_time (seconds)

    Used for waiting for qsub to finish
    """
    print("GOT HERE")
    time_in_loop = 0
    while not os.path.exists(file_path):
        if time_in_loop < wait_time:
            print("waiting")
            time.sleep(1)
            time_in_loop += 1
        else:
            raise IOError(
                "Cannot find file {} within {} seconds".format(file_path, wait_time)
            )


def chunks(l, n):
    """ Divide a list l into chunks of length n.

    Yield with consecutive letters

    Used for splitting a list of atomic points into breaks of 9999 for
    display in a pdb file.

    Parameters
    -------------------
    l: list
        list to be split into chunks
    y: int
        length to chunk list into

    Yields
    ----------------------
    list
        slice of original list up to chunk size
    str
        A letter associated with chunk (A-Z)
    """

    alphabet = []
    for letter in range(65, 91):
        alphabet.append(chr(letter))

    # For item i in a range that is a length of l,
    pos = 0
    for i in xrange(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i : i + n], alphabet[pos]
        pos += 1


def append_csv(in_csv1, in_csv2, out_csv):
    """
    Append two input csv files to out_csv

    Parameters
    ----------
    in_csv1: str
        path to input csv 1
    in_csv2: str
        path to input csv 2
    out_csv: str
        path to output csv

    Returns
    -------
    None

    Notes
    -----
    Writes csv to out_csv
    """

    df1 = pd.read_csv(in_csv1, header=None)
    df2 = pd.read_csv(in_csv2, header=None)

    df = df1.append(df2, ignore_index=True)

    df.to_csv(out_csv, header=None, index=False)
