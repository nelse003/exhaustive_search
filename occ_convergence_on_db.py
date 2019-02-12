import pandas as pd
import numpy as np
import time
import iotbx.pdb
import sys
import traceback
import itertools

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

def drop_similiar_cols(cols, df, atol=0.001):

    print("COLS:{}".format(cols))
    print(df)

    for pair_col in itertools.combinations(cols, 2):

        if len(pair_col) == 2:
            col_1 = df[pair_col[0]]
            col_2 = df[pair_col[1]]

            if np.allclose(col_1.values, col_2.values, atol=atol):
                cols_to_drop = pair_col[0]
                return df.drop(cols_to_drop, axis=1)
            else:
                return df
        else:
            return df

def get_occupancy_df(log_path, pdb_path, crystal, lig_name='LIG',):

    chain = get_lig_chain_id(pdb_path=pdb_path, lig_name=lig_name)

    lig_groups = read_lig_occ_group_from_refmac_log(log_path=log_path,
                                                    lig_chain=chain)
    lig_groups = np.asarray(lig_groups, dtype=np.int64)

    occ_conv_df = read_occupancies_from_refmac_log(log_path)

    # Get columns associated with bound & ground from
    # ligand groups

    column_bound = list(
        set(occ_conv_df.columns.values).intersection(
            set(lig_groups)))

    if len(column_bound) == 0:
        raise ValueError("No bound state identified")

    column_ground = list(
        set(occ_conv_df.columns.values).symmetric_difference(
            set(lig_groups)))

    if len(column_bound) == 0:
        raise ValueError("No ground state identified")

    # Drop duplicate columns
    short_occ_conv_df = drop_similiar_cols(column_ground,
                                        occ_conv_df,
                                        atol=0.001)
    short_occ_conv_df = drop_similiar_cols(column_bound,
                                           short_occ_conv_df,
                                        atol=0.001)

    # Return occupancy to full
    multiplicity_of_occ_groups = len(occ_conv_df.T) / \
                                 len(short_occ_conv_df.T)
    occ_corrected_df = short_occ_conv_df * multiplicity_of_occ_groups

    # Check whether all occupancy values are 1.0
    if len(pd.unique(occ_corrected_df.values)) == 1:
        if pd.unique(occ_corrected_df.values)[0] == 1.0:
            raise ValueError("All occupancies are {}. "
                             "Likely there is no ground state".format(
                             pd.unique(occ_corrected_df.values)[0]))
        else:
            print(occ_corrected_df)
            print(occ_conv_df)
            pass

    if len(occ_corrected_df.T) != 2:
        print(occ_conv_df)
        print(occ_corrected_df)
        raise ValueError("Too many columns in occupancy df: {}".format(occ_corrected_df))

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

    occ_df = occ_corrected_df.rename(columns={column_bound: "bound",
                                              column_ground: "ground"},
                                              index=str)
    # Transposing and adding crystal_id
    occ_df = occ_df.T
    occ_df.reset_index(level=0, inplace=True)
    occ_df.rename(columns={'index': 'state'},
                  index=str,
                  inplace=True)
    occ_df['crystal'] = [crystal, crystal]
    occ_df = occ_df.set_index(['crystal', 'state'])

    return occ_df

if __name__ == "__main__":
    """
    Process log files from refmac runs to get occupancy convergence
    
    Requires ccp4-python
    
    """

    # TODO Sort utils/plotting into ccp4 dependent and ccp4 non dependent sections

    log_csv = "/dls/science/groups/i04-1/elliot-dev/Work/" \
              "exhaustive_parse_xchem_db/log_pdb_mtz.csv"

    occ_conv_csv = "/dls/science/groups/i04-1/elliot-dev/Work/" \
              "exhaustive_parse_xchem_db/occ_conv.csv"

    occ_conv_fails_csv = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_parse_xchem_db/occ_conv_failures.csv"

    log_df = pd.read_csv(log_csv)

    #individual test case where columns
    # should overlap but don't because
    # non-identical
    #
    # merge_log = "/dls/labxchem/data/2018/lb19005-1/processing/" \
    # "analysis/initial_model/PaFEN_C3-x0440/" \
    # "Refine_0001/refine_1.quick-refine.log"
    #
    # merge_pdb = "/dls/labxchem/data/2018/lb19005-1/processing/" \
    #             "analysis/initial_model/PaFEN_C3-x0440/" \
    #             "Refine_0001/refine_1.pdb"
    #
    # crystal, target = pdb_path_to_crystal_target(merge_pdb)
    # occ_df = get_occupancy_df(log_path=merge_log,
    #                           pdb_path=merge_pdb,
    #                           lig_name='LIG',
    #                           crystal=crystal)


    occ_conv_df_list = []
    failures = []
    for index, row in log_df.iterrows():
        print(row.refine_log)

        try:
            crystal, target = pdb_path_to_crystal_target(row.pdb_latest)
            occ_df = get_occupancy_df(log_path=row.refine_log,
                                      pdb_path=row.pdb_latest,
                                      lig_name='LIG',
                                      crystal=crystal)
        except ValueError:
            # This is for handling the exception,
            # utilising it's traceback to determine
            # the consitent errors
            #
            # This is python 2.7 specific
            #
            # Could alternatively be done by logging

            ex_type, ex, tb = sys.exc_info()
            tb_txt = traceback.extract_tb(tb)
            error = (row.refine_log, ex_type, ex, tb_txt)
            failures.append(error)
            print('Error')
            continue

        occ_conv_df_list.append(occ_df)

        # os.makedirs()
        #
        # plot_occupancy_convergence(occ_conv_df=occ_df,
        #                            plot_filename=
        #                            "/dls/science/groups/i04-1/elliot-dev/"
        #                            "Work/exhaustive_parse_xchem_db/"
        #                            "log_png_test.png")

    failures_df =pd.DataFrame(failures,columns=['log',
                                                'exception_type',
                                                'excpetion',
                                                'traceback'])
    failures_df.to_csv(occ_conv_fails_csv)

    occ_conv_summary_df = pd.concat(occ_conv_df_list)
    occ_conv_summary_df.to_csv(occ_conv_csv)