import os
import pandas as pd

from exhaustive.exhaustive import master_phil


# from giant.jiffies.score_model import run as score_model
# from giant.jiffies.score_model import master_phil as score_phil
# import libtbx.phil

# Changing project just requires changing the directory.
# Note the /mnt prefix for running locally on laptop

# ini_folder = "/mnt/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios"
# ini_folder = "/mnt/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/titration_series"

# This works in ccp4 python. Currently using the score model script as the
# score params stuff is only needed if changing to a different residue
# such as GDP for janine


def get_xtals_from_db(
    params,
    refinement_outcomes="'3 - In Refinement',"
    "'4 - CompChem ready', "
    "'5 - Deposition ready',"
    "'6 - Deposited'",
):
    assert os.path.isfile(
        params.input.database_path
    ), "The database file: \n {} \n does not exist".format(params.input.database_path)

    # Open connection to sqlite database
    conn = sqlite3.connect(params.input.database_path)
    cur = conn.cursor()

    cur.execute(
        "SELECT CrystalName, RefinementPDB_latest, RefinementMTZ_latest "
        "FROM mainTable WHERE RefinementOutcome in ({})"
        " AND  (RefinementPDB_latest AND RefinementMTZ_latest) IS NOT NULL".format(
            refinement_outcomes
        )
    )

    refinement_xtals = cur.fetchall()

    # Close connection to the database
    cur.close()

    for xtal_name, pdb, mtz in refinement_xtals:
        pdb = pdb.encode("ascii")
        mtz = mtz.encode("ascii")
        xtal_name = xtal_name.encode("ascii")
        yield xtal_name, pdb, mtz


# params.input.in_path = "/dls/labxchem/data/2016/lb13385-64/processing/analysis/initial_model"
params = master_phil.extract()
#
params.input.database_path = (
    "/dls/labxchem/data/2016/lb13385-64/processing/database/soakDBDataFile.sqlite"
)
params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus"

for xtal_name, pdb, mtz in get_xtals_from_db(
    params,
    refinement_outcomes="'4 - CompChem ready', "
    "'5 - Deposition ready',"
    "'6 - Deposited'",
):

    os.system(
        "cp {} {}".format(
            pdb, os.path.join(params.output.out_dir, xtal_name, os.path.basename(pdb))
        )
    )
    os.system(
        "cp {} {}".format(
            mtz, os.path.join(params.output.out_dir, xtal_name, os.path.basename(mtz))
        )
    )
    os.chdir(os.path.join(params.output.out_dir, xtal_name))

    if not os.path.exists(
        os.path.join(params.output.out_dir, xtal_name, "residue_scores.csv")
    ):
        os.system("giant.score_model {} {}".format(pdb, mtz))

# Currently pandas is failing to import in ccp4 python so this is done seperately in a conda env

dfs = []

for xtal_name, pdb, mtz in get_xtals_from_db(
    params,
    refinement_outcomes="'4 - CompChem ready', "
    "'5 - Deposition ready',"
    "'6 - Deposited'",
):
    edstats_csv = os.path.join(params.output.out_dir, xtal_name, "residue_scores.csv")
    edstats_df = pd.read_csv(edstats_csv)
    edstats_df["Dataset"] = xtal_name
    dfs.append(edstats_df)

print("-------------------------------------------")
compound_edstats = pd.concat(dfs, ignore_index=True)
print(compound_edstats)
compound_edstats.to_csv(os.path.join(params.output.out_dir, "edstats.csv"))
