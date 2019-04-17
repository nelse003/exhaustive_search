import sys
import os
import time
import pandas as pd

from exhaustive import master_phil
from exhaustive import run as exhaustive
from utils.utils import u_iso_to_b_fac, get_minimum_fofc


def append_csv(in_csv1, in_csv2, out_csv):
    """
    Append two csv to out_csv

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
    """

    df1 = pd.read_csv(in_csv1, header=None)
    df2 = pd.read_csv(in_csv2, header=None)

    df = df1.append(df2, ignore_index = True)

    df.to_csv(out_csv, header=None, index=False)


def run(params):
    """
    Run Exhaustive search twice at 0.05 and 0.01 spacing

    Run over whole occupancy range [0,1] at 0.05 spacing,
    then rerun near minimum at 0.01 spacing [min - 0.1, min + 0.1].

    Parameters
    ----------
    params

    Returns
    -------
    None

    Notes
    -----
    Writes csv to params.exhaustive.output.csv_name
    """
    base_csv = params.exhaustive.output.csv_name.split('.csv')[0]

    params.exhaustive.options.step = 0.05
    params.exhaustive.output.csv_name = base_csv + "_coarse.csv"

    t1 = time.time()
    exhaustive(params)
    t2 = time.time()

    occ, u_iso, fo_fc = get_minimum_fofc(
        os.path.join(params.output.out_dir,
                     params.exhaustive.output.csv_name))

    occ_min = occ - 0.1
    occ_max = occ + 0.1
    u_iso_min = u_iso - 0.1
    u_iso_max = u_iso + 0.1

    params.exhaustive.options.step = 0.01

    params.exhaustive.output.csv_name = base_csv + "_fine.csv"

    params.exhaustive.options.lower_occ = occ_min
    params.exhaustive.options.upper_occ = occ_max
    params.exhaustive.options.lower_u_iso = u_iso_min
    params.exhaustive.options.upper_u_iso = u_iso_max

    t3 = time.time()
    exhaustive(params)
    t4= time.time()

    occ_fine, u_iso_fine, fo_fc_fine = get_minimum_fofc(os.path.join(params.output.out_dir,
                                                     params.exhaustive.output.csv_name))

    print("Minimum from first run is occupancy: {occupancy},\n"
          "B factor {b_factor}\n"
          "This took {time} seconds".format(occupancy=occ,
                                            b_factor=u_iso_to_b_fac(u_iso),
                                            time=t2-t1))

    print("Minimum from second run is occupancy: {occupancy},\n"
          "B factor {b_factor}\n"
          "This took {time} seconds".format(occupancy=occ_fine,
                                         b_factor=u_iso_to_b_fac(u_iso_fine),
                                         time=t4-t3))

    append_csv(in_csv1=os.path.join(params.output.out_dir,
                                    base_csv + "_coarse.csv"),
               in_csv2=os.path.join(params.output.out_dir,
                                    base_csv + "_fine.csv"),
               out_csv=os.path.join(params.output.out_dir,
                                    base_csv +".csv"))

if(__name__ == "__main__"):
    from giant.jiffies import run_default

    blank_arg_prepend = {'.pdb': 'pdb=', '.mtz': 'mtz=', '.csv': 'csv='}

    run_default(run=run, master_phil=master_phil,
                blank_arg_prepend=blank_arg_prepend, args=sys.argv[1:])
