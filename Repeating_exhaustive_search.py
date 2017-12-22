import os
import sys
import libtbx.phil
import pandas as pd
import csv
from exhaustive_search import cmd_run as exhaustive_search, get_minimum_fofc
from plot_exhaustive_search import scatter_plot_4col

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

def run(params):

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
    with open("min_occ_u_iso",'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        for xtal_name, pdb, mtz in parse_repeat_soak_csv(params):
            if pdb and mtz is not None:
                try:
                    assert os.path.exists(pdb), 'PDB File does not exist: {}'.format(pdb)
                    assert os.path.exists(mtz), 'MTZ File does not exist: {}'.format(mtz)
                    os.chdir("output_repeat_soak/{}".format(xtal_name))
                    occ, u_iso = get_minimum_fofc("{}_0A_bound_ground_covary_frac_fix".format(xtal_name))
                    row = [xtal_name, occ, u_iso]
                    writer.writerow(row)
                    sys.stdout.flush()
                    os.chdir("../..")
                except:
                    continue
if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)


