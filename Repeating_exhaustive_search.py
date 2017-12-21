import os
import libtbx.phil
import pandas as pd

##############################################################

PROGRAM = 'Exhaustive Search'
DESCRIPTION = """
    Take in pdb and mtz, or csv of pdb and mtz locations,
    run a exhaustive search of B factor and occupancy to determine unique minima.
"""
blank_arg_prepend = {'.pdb': 'pdb=','.mtz'}

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
        yield row["RefinementPDB_latest"], row["RefinementMTZ_latest"]

def run(params):

    if not os.path.exists(params.input.csv):
        assert os.path.exists(params.input.pdb), 'PDB File does not exist'
        assert os.path.exists(params.input.mtz), 'MTZ File does not exist'
        os.system("exhaustive_search.py {} {}".format(pdb,mtz))

    elif os.path.exists(params.input.csv):
        for pdb, mtz in parse_repeat_soak_csv(params):
            os.system("exhasutive_search.py {} {}".format(pdb,mtz))

    else:
        print ("Please supply a pdb and mtz, or a csv file")

if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)


