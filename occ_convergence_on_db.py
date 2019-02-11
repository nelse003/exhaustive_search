import pandas as pd
import numpy as np
import time
import iotbx.pdb

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

def read_lig_occ_group_from_refmac_log(log_path, lig_chain):

    lig_occ_groups = []
    with open(log_path, 'r') as log:
        for line in log:
            if "occupancy group" in line:
                if "chain {}".format(lig_chain) in line:
                    if "Data line" not in line:
                        lig_occ_groups.append(line.split('id')[1].split()[0])

    return lig_occ_groups


def pdb_path_to_crystal_target(pdb_path):

    pdb_path_split = pdb_path.split('/')
    crystal_name = [path_part for path_part in pdb_path_split if '-x' in path_part]
    crystal_name = crystal_name[0]
    target = crystal_name.split('-')[0]

    return crystal_name, target

def get_lig_chain_id(pdb_path, lig_name='LIG'):
    """Get chain which ligand is in
    
    Notes
    --------------
    Requires ccp4-python (iotbx.pdb)
    """

    pdb_in = iotbx.pdb.hierarchy.input(file_name=pdb_path)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()
    lig_sel = sel_cache.selection("resname {}".format(lig_name))
    lig_hierarchy = pdb_in.hierarchy.select(lig_sel)
    lig_chains = []
    for chain in lig_hierarchy.only_model().chains():
        lig_chains.append(chain.id)

    if len(lig_chains) == 1:
        return chain.id
    else:
        raise ValueError('Ligand chain is not unique: {}'.format(lig_chains))

def get_occupancy_df(log_path, pdb_path, lig_name='LIG'):

    chain = get_lig_chain_id(pdb_path=pdb_path, lig_name=lig_name)

    lig_groups = read_lig_occ_group_from_refmac_log(log_path=log_path,
                                                    lig_chain=chain)
    lig_groups = np.asarray(lig_groups, dtype=np.int64)

    occ_conv_df = read_occupancies_from_refmac_log(log_path)

    # Get columns associated with bound
    column_bound = list(
        set(occ_conv_df.columns.values).intersection(
            set(lig_groups)))

    column_ground = list(
        set(occ_conv_df.columns.values).symmetric_difference(
            set(lig_groups)))

    # Remove duplicate columns
    short_occ_conv_df = occ_conv_df.T.drop_duplicates().T
    multiplicity_of_occ_groups = len(occ_conv_df.T) / \
                                 len(short_occ_conv_df.T)
    occ_corrected_df = short_occ_conv_df * multiplicity_of_occ_groups

    if len(occ_corrected_df.T) != 2:
        raise ValueError("Too many columns in occupancy df: "
                         "{}".format(occ_corrected_df))

    # Handle labels for duplicate columns
    for col in column_ground:
        if col not in occ_corrected_df.columns.values:
            column_ground.remove(col)

    if len(column_ground) != 1:
        raise ValueError("Unexpected number of "
                         "ground-state associated "
                         "columns: {}".format(column_ground))
    column_ground = np.int64(column_ground[0])

    for col in column_bound:
        if col not in occ_corrected_df.columns.values:
            column_bound.remove(col)

    if len(column_bound) != 1:
        raise ValueError("Unexpected number of "
                         "bound-state associated "
                         "columns: {}".format(column_bound))
    column_bound = np.int64(column_bound[0])

    occ_corrected_df = occ_corrected_df.rename(columns={column_bound: "bound",
                                                        column_ground: "ground"},
                                               index=str)

    return occ_corrected_df

if __name__ == "__main__":
    """
    Process log files from refmac runs to get occupancy convergence
    
    #TODO Sort utils/plotting into ccp4 dependent 
    and ccp4 non dependent sections
    
    Currently requires ccp4-python
    
    """


    log_csv = "/dls/science/groups/i04-1/elliot-dev/Work/" \
              "exhaustive_parse_xchem_db/log_pdb_mtz.csv"

    occ_conv_csv = "/dls/science/groups/i04-1/elliot-dev/Work/" \
              "exhaustive_parse_xchem_db/occ_conv.csv"

    log_df = pd.read_csv(log_csv)

    occ_conv_df_list = []
    failures = []
    for index, row in log_df.iterrows():


        print(row.refine_log)
        try:
            crystal, target = pdb_path_to_crystal_target(row.pdb_latest)
            occ_df = get_occupancy_df(log_path=row.refine_log,
                                      pdb_path=row.pdb_latest,
                                      lig_name='LIG')
        # TODO handle exceptions saving tracebakcs?
        except (ValueError, AssertionError, IndexError) as e:
            failures.append(row.refine_log)
            continue

        occ_conv_df_list.append(occ_df)

        # os.makedirs()
        #
        # plot_occupancy_convergence(occ_conv_df=occ_df,
        #                            plot_filename=
        #                            "/dls/science/groups/i04-1/elliot-dev/"
        #                            "Work/exhaustive_parse_xchem_db/"
        #                            "log_png_test.png")

    # TODO Transpose and add xtal_id?
    occ_conv_summary_df = pd.concat(occ_conv_df_list)
    occ_conv_summary_df.to_csv(occ_conv_csv)
    print(failures)
    print(len(failures))