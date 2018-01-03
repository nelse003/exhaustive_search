from __future__ import print_function
import os
import sys
import libtbx.phil
import pandas as pd
import csv
from exhaustive_search import cmd_run as exhaustive_search, get_minimum_fofc
from plot_exhaustive_search import scatter_plot_4col, scatter_plot
import sqlite3

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
}
output{
    out_dir = "/hdlocal/home/enelson/Dropbox/DPhil/exhaustive_search/output/"
        .type = str
}
options{

}
""", process_includes=True)

########################################################################
def parse_repeat_soak_csv(params):

    input_df = pd.read_csv(params.input.csv)
    for index, row in input_df.iterrows():
        yield row["CrystalName"],row["RefinementPDB_latest"], row["RefinementMTZ_latest"]


def get_in_refinement_or_better(params):

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

    # TODO Move repat soaks to function
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

    # TODO Move saving all minima to function
    # Saving all minima

    # with open("min_occ_u_iso",'w') as f1:
    #     writer = csv.writer(f1, delimiter=',', lineterminator='\n')
    #     for xtal_name, pdb, mtz in parse_repeat_soak_csv(params):
    #         if pdb and mtz is not None:
    #             try:
    #                 assert os.path.exists(pdb), 'PDB File does not exist: {}'.format(pdb)
    #                 assert os.path.exists(mtz), 'MTZ File does not exist: {}'.format(mtz)
    #                 os.chdir("output_repeat_soak/{}".format(xtal_name))
    #                 occ, u_iso = get_minimum_fofc("{}_0A_bound_ground_covary_frac_fix".format(xtal_name))
    #                 row = [xtal_name, occ, u_iso]
    #                 writer.writerow(row)
    #                 sys.stdout.flush()
    #                 os.chdir("../..")
    #             except:
    #                 continue

    for xtal_name, pdb, mtz in get_in_refinement_or_better(params):


        assert os.path.exists(pdb), 'PDB File does not exist: {}'.format(pdb)
        assert os.path.exists(mtz), 'MTZ File does not exist: {}'.format(mtz)

        #### For Exhaustive search run ####

        # args = [pdb, mtz]
        # exhaustive_search(args, xtal_name)

        #### For Plotting ####

        # output_folder = "output_DCP2_refinements/{}".format(xtal_name)
        # output_path = os.path.join(os.getcwd(), output_folder)
        # output_path_base = os.path.join(os.getcwd(), "output_DCP2_refinements")
        #
        # if not os.path.exists(output_path_base):
        #     os.mkdir(output_path_base)
        #
        # if not os.path.exists(output_path):
        #     os.mkdir(output_path)
        # os.chdir(output_folder)
        # csv_name = '{}_0A_bound_ground_covary_frac_fix'.format(xtal_name)
        # scatter_plot_4col(csv_name)
        # os.chdir("../../")

        #### For getting single plot #####
        if xtal_name == "DCP2B-x0020":
            #args = [pdb, mtz]
            #exhaustive_search(args, xtal_name)
            csv_name = 'mean_covary_0A_buffer'
            output_folder = "output_DCP2_refinements/{}".format(xtal_name)
            output_path = os.path.join(os.getcwd(), output_folder)
            output_path_base = os.path.join(os.getcwd(), "output_DCP2_refinements")

            if not os.path.exists(output_path_base):
                os.mkdir(output_path_base)

            if not os.path.exists(output_path):
                os.mkdir(output_path)
            print("AAAAAAAAAAA")
            os.chdir(output_folder)
            print(os.getcwd())
            scatter_plot(csv_name)
            sys.exit()







if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)


