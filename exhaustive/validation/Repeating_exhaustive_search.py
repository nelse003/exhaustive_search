from __future__ import print_function

import os
import sqlite3
import sys

import libtbx.phil
import pandas as pd

from exhaustive.exhaustive import run as exhaustive_search

##############################################################

PROGRAM = 'Repeat Exhaustive Search'
DESCRIPTION = """
    Take in csv of pdb and mtz locations,
    run a exhaustive search of B factor and occupancy to determine unique minima.
"""
blank_arg_prepend = {'.pdb': 'pdb=','.mtz': 'mtz=','.csv': 'csv='}

##############################################################
master_phil = libtbx.phil.parse("""
input{
    pdb = None
        .type = path
    mtz = None
        .type = path
    csv = None
        .type = path
    database_path = None
        .type = path
    csv_name = 'u_iso_occupancy_vary_new_atoms'
        .type = str
}
output{
    out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/occupancy_group_with_refinement"
        .type = str
    minima_csv_name = "min_occ_u_iso_all"
        .type = str
}
options{

}
""", process_includes=True)
########################################################################
import logging
import datetime

logging.basicConfig(filename=datetime.datetime.now().strftime('/dls/science/groups/i04-1/elliot-dev/Work/' \
                                                              'exhaustive_search/logs/exhaustive_search_%H_%M_%d_%m_%Y.log'),
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


# # This is used to allow any excpetion, such as a failed assert statement to be logged, doesn't also go to sysout.
# def excepthook(*args):
#   logging.getLogger().error('Uncaught exception:', exc_info=args)
# sys.excepthook = excepthook

########################################################################
def parse_repeat_soak_csv(params):

    input_df = pd.read_csv(params.input.csv)
    for index, row in input_df.iterrows():
        yield row["CrystalName"],row["RefinementPDB_latest"], row["RefinementMTZ_latest"]

def get_in_refinement_or_better(params):
    assert os.path.isfile(params.input.database_path), \
        "The database file: \n {} \n does not exist".format(params.input.database_path)

    # Open connection to sqlite database
    conn = sqlite3.connect(params.input.database_path)
    cur = conn.cursor()

    cur.execute("SELECT CrystalName, RefinementPDB_latest, RefinementMTZ_latest "
                "FROM mainTable WHERE RefinementOutcome in ('3 - In Refinement',"
                "'4 - CompChem ready','5 - Deposition ready','6 - Deposited')" 
                " AND  (RefinementPDB_latest AND RefinementMTZ_latest) IS NOT NULL")

    refinement_xtals = cur.fetchall()

    # Close connection to the database
    cur.close()

    for xtal_name, pdb, mtz in refinement_xtals:
        pdb = pdb.encode('ascii')
        mtz = mtz.encode('ascii')
        xtal_name = xtal_name.encode('ascii')
        yield xtal_name, pdb, mtz


def run(params):

    #get_all_minima(params)

    # Repeat soaks of DCP2B; run over all

    # if not os.path.exists(params.input.csv):
    #     assert os.path.exists(params.input.pdb), 'PDB File does not exist: {}'.foramt(params.input.pdb)
    #     assert os.path.exists(params.input.mtz), 'MTZ File does not exist: {}'.format(params.input.mtz)
    #     exhaustive_search(args, xtal_name)
    #
    # elif os.path.exists(params.input.csv):
    #     for xtal_name, pdb, mtz in parse_repeat_soak_csv(params):
    #         if pdb and mtz is not None:
    #             try:
    #                 assert os.path.exists(pdb), 'PDB File does not exist: {}'.format(pdb)
    #                 assert os.path.exists(mtz), 'MTZ File does not exist: {}'.format(mtz)
    #                 args = [pdb, mtz]
    #                 exhaustive_search(args, xtal_name)
    #             except:
    #                 print "Skipping"
    #                 continue
    #         else:
    #             print "No pdb/mtz combo for this repeat, contuining"
    #             continue
    # else:
    #     print ("Please supply a pdb and mtz, or a csv file")


    if not os.path.exists(params.output.out_dir):
        logger.info('Creating output directory {}'.format(params.output.out_dir))
        os.mkdir(params.output.out_dir)
    else:
        logger.info('Output directory {} exists and is being used'.format(params.output.out_dir))

    logger.info('Looping over all files that are \'in refinement\' '
                'or better in the supplied datafile: \n {}'.format(params.input.database_path))



    for xtal_name, pdb, mtz in get_in_refinement_or_better(params):

        logger.info(xtal_name)

        assert os.path.exists(pdb), 'PDB File does not exist: {}'.format(pdb)
        assert os.path.exists(mtz), 'MTZ File does not exist: {}'.format(mtz)

        os.chdir(os.path.join(params.output.out_dir))

        #### For Exhaustive search run ####
        args = [pdb, mtz]

        try:
            exhaustive_search(args, xtal_name)
        except UnboundLocalError:
            logger.info("Skipping onto the next crystal")
            continue

        if not os.path.exists(os.path.join(params.output.out_dir, xtal_name)):
            os.mkdir(os.path.join(params.output.out_dir, xtal_name))
            os.chdir(os.path.join(params.output.out_dir, xtal_name))
        else:
            os.chdir(os.path.join(params.output.out_dir, xtal_name))
        #scatter_plot(params.input.csv_name)

        logger.info('Completed: {}'.format(xtal_name))
        #



    #
    #     #### For Plotting ####
    #
    #     output_folder = "test_runs/{}".format(xtal_name)
    #     output_path = os.path.join(os.getcwd(), output_folder)
    #     output_path_base = os.path.join(os.getcwd(), "output_DCP2_refinements")
    #
    #     if not os.path.exists(output_path_base):
    #         os.mkdir(output_path_base)
    #
    #     if not os.path.exists(output_path):
    #         os.mkdir(output_path)
    #     os.chdir(output_folder)
    #     csv_name = 'u_iso_occupancy_vary'
    #     print(os.getcwd())
    #     os.rename(csv_name, csv_name + ".csv")
    #     scatter_plot(csv_name)
            #     os.chdir("../../"


if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)


