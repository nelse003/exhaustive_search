import argparse
import os
from shutil import rmtree
import pandas as pd
import sqlalchemy

pd.options.display.max_columns = 30
pd.options.display.max_colwidth = 50

def parse_args():
    parser = argparse.ArgumentParser()

    # PanDDA run
    parser.add_argument("-r", "--pandda_run")
    # default="/dls/labxchem/data/2018/lb19758-9/processing/analysis/panddas_run12_all_onlydmso"

    # Output
    parser.add_argument("-o", "--output", default="/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_parse_xchem_db")

    # Database
    parser.add_argument("-hs", "--host", default="172.23.142.43")
    parser.add_argument("-p", "--port", default="5432")
    parser.add_argument("-d", "--database", default="test_xchem")
    parser.add_argument("-u", "--user", default="conor")
    parser.add_argument("-pw", "--password", default="c0n0r")

    # Interpeter args
    #parser.add_argument("-cpy", "--ccp4_python", default="/home/zoh22914/myccp4python")

    # Operation
    parser.add_argument("-b", "--build", default=0)


    args = parser.parse_args()

    return args

def get_databases(args):

    engine=sqlalchemy.create_engine("postgresql://{}:{}@{}:{}/{}".format(
                                    args.user, args.password, args.host, 
                                    args.port, args.database))

    print(engine.table_names())

    pandda_analysis = pd.read_sql_query("SELECT * FROM pandda_analysis", 
                                        con=engine)

    pandda_run = pd.read_sql_query("SELECT * FROM pandda_run", con=engine)

    pandda_events = pd.read_sql_query("SELECT * FROM pandda_event", con=engine)

    pandda_events_stats = pd.read_sql_query("SELECT * FROM pandda_event_stats", con=engine)

    statistical_maps = pd.read_sql_query("SELECT * FROM pandda_statistical_map", con=engine)

    crystals = pd.read_sql_query("SELECT * FROM crystal", con=engine)

    ligands = pd.read_sql_query("SELECT * FROM proasis_hits", con=engine)

    ligand_stats = pd.read_sql_query("SELECT * FROM ligand_edstats", con=engine)

    data_processing = pd.read_sql_query("SELECT * FROM data_processing", con=engine)

    refinement = pd.read_sql_query("SELECT * FROM refinement",con=engine)

    compounds = pd.read_sql_query("SELECT * FROM compounds", con=engine)

    databases = {"pandda_analyses": pandda_analysis,
                 "pandda_runs": pandda_run ,
                 "pandda_events": pandda_events,
                 "pandda_event_stats": pandda_events_stats,
                 "statistical_maps": statistical_maps,
                 "crystals": crystals,
                 "ligands": ligands,
                 "ligand_stats": ligand_stats,
                 "data_processing": data_processing,
                 "refinement": refinement,
                 "compounds": compounds}

    return databases

if __name__ == "__main__":

    # get args
    print("Parsing arguments...")
    args = parse_args()

    # Set up output dir
    print("Establishing output directory...")
    try:
        os.mkdir(args.output)
    except:
        rmtree(args.output)
        os.mkdir(args.output)

    # get_databases
    print("Setting up databases...")
    databases = get_databases(args)

    #print(databases['refinement'].head())

    refine_df = databases['refinement']

    # Refinemnt columns
    #
    # ['id' 'bound_conf' 'cif' 'cif_prog' 'cif_status' 'lig_bound_conf' 'lig_cc'
    # 'lig_confidence' 'matrix_weight' 'molprobity_score' 'mtz_free'
    # 'mtz_latest' 'outcome' 'pdb_latest' 'r_free' 'ramachandran_favoured'
    # 'ramachandran_outliers' 'rcryst' 'refinement_path' 'res' 'rmsd_angles'
    # 'rmsd_bonds' 'spacegroup' 'status' 'crystal_name_id' 'lig_confidence_int'
    # 'lig_confidence_string']

    # TODO Check columns where bound conf does not exist to whether refine.pdb exists

    print("Length of refinement table: {}".format(len(refine_df)))

    # Check rows with just superposed confirmation
    superposed_df = refine_df[refine_df.pdb_latest.notnull()]

    print("Length of refinement table with superposed state"
          "specified: {}".format(len(superposed_df)))

    # Get only rows with bound conformation
    bound_df = refine_df[refine_df.bound_conf.notnull()]

    print("Length of refinement table with bound state"
          "specified: {}".format(len(bound_df)))

    # drops NaN if both bound_conf and pdb_latest are NaN
    with_pdb_df = refine_df.dropna(axis='index',how='all', subset=['bound_conf','pdb_latest'])

    print("Length of refinement table with either bound state"
          "or superposed state specified: {}".format(len(with_pdb_df)))

    # Check path exists for bound, ground and refine.pdb
    bound_missing_ids = []
    superposed_missing_ids = []
    ground_missing_ids = []

    for index, row in with_pdb_df.iterrows():

        if row.bound_conf is not None:
            bound_conf_dir = os.path.dirname(row.bound_conf)
            bound_conf_file_name = os.path.basename(row.bound_conf)
            ground_conf_file_name = bound_conf_file_name.replace('bound','ground')
            ground_conf = os.path.join(bound_conf_dir,ground_conf_file_name)

            if not os.path.isfile(row.bound_conf):
                bound_missing_ids.append(row.id)

            if not os.path.isfile(ground_conf):
                ground_missing_ids.append(row.id)
        else:
            bound_missing_ids.append(row.id)
            ground_missing_ids.append(row.id)

        if row.pdb_latest is not None:
            if not os.path.isfile(row.pdb_latest):
                superposed_missing_ids.append(row.id)
        else:
            superposed_missing_ids.append(row.id)


    print("Number of rows where bound structure is not present in filesystem: "
          "{}".format(len(bound_missing_ids)))

    print("Number of rows where superposed structure is not present in "
          "filesystem: {}".format(len(superposed_missing_ids)))

    print("Number of rows where ground structure is not present in "
          "filesystem: {}".format(len(ground_missing_ids)))



    # print(len(bound_conf_df))
    # print(refine_df.columns.values)
    #


