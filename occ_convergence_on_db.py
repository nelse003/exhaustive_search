import pandas as pd
import numpy as np
import time
from exhaustive.exhaustive.plotting.plot import plot_occupancy_convergence

def read_occupancies_from_refmac_log(log_path):
    """Read group occupancies from log

    Line to be read from is Group occupances  0.12964654  0.122 ...

    Parameters
    -----------
    log_path: str
        Path to log file

    Returns
    -----------
    A dataframe of columns = cycles, rows= occ groups

    """
    group_occupancies = {}
    with open(log_path, 'r') as log:
        for line in log:
            if 'CGMAT cycle number =' in line:
                cycle_line = line
                cycle_int = int(cycle_line.split('=')[-1].lstrip(' ').rstrip('\n'))
            if 'Group occupances' in line:
                group_occupancies[cycle_int] = [float(group_occ) for group_occ in line.split()[2:]]

    df = pd.DataFrame.from_dict(group_occupancies, orient='index')
    df.columns = np.arange(1, len(df.columns) + 1)

    return df

if __name__ == "__main__":

    log_csv = "/dls/science/groups/i04-1/elliot-dev/Work/" \
              "exhaustive_parse_xchem_db/log_pdb_mtz.csv"

    log_df = pd.read_csv(log_csv)

    for index, row in log_df.iterrows():

        start = time.time()
        df = read_occupancies_from_refmac_log(row.refine_log)
        end = time.time()
        print(end - start)
        start = time.time()
        plot_occupancy_convergence(occ_conv_df=df,
                                   plot_filename="/dls/science/groups/i04-1/elliot-dev/Work/" \
              "exhaustive_parse_xchem_db/log_png_test.png")
        end = time.time()
        print(end - start)
        break