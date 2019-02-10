import argparse
import os
from shutil import rmtree
import pandas as pd
import sqlalchemy


def parse_args():
    parser = argparse.ArgumentParser()

    # PanDDA run
    parser.add_argument("-r", "--pandda_run")
    # default="/dls/labxchem/data/2018/lb19758-9/processing/analysis/panddas_run12_all_onlydmso"

    # Output
    parser.add_argument("-o", "--output", default="harness_test/")

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

    engine=sqlalchemy.create_engine("postgresql://{}:{}@{}:{}/{}".format(args.user, args.password, args.host, args.port, args.database))

    pandda_analysis = pd.read_sql_query("SELECT * FROM pandda_analysis", con=engine)

    pandda_run = pd.read_sql_query("SELECT * FROM pandda_run", con=engine)

    pandda_events = pd.read_sql_query("SELECT * FROM pandda_event", con=engine)

    pandda_events_stats = pd.read_sql_query("SELECT * FROM pandda_event_stats", con=engine)

    statistical_maps = pd.read_sql_query("SELECT * FROM pandda_statistical_map", con=engine)

    crystals = pd.read_sql_query("SELECT * FROM crystal", con=engine)

    ligands = pd.read_sql_query("SELECT * FROM proasis_hits", con=engine)

    ligand_stats = pd.read_sql_query("SELECT * FROM ligand_edstats", con=engine)

    data_processing = pd.read_sql_query("SELECT * FROM data_processing", con=engine)

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
                 "compounds": compounds}

    return databases


if __name__ == "__main__":

    # get args
    print("Parsing arguments...")
    args = parse_args()

    # Constants
    # print("Defining constants...")
    # home = "/home/zoh22914/"
    # base_python = "python"
    # ccp4_python = home + "myccp4python"
    # np2ccp4 = home + "pandda-conor/lib-python/ppandda/investigations/numpy_to_ccp4.py"
    # ccp42np = home + "pandda-conor/lib-python/ppandda/investigations/ccp4_to_np.py"

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