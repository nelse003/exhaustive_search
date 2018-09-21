import os
import pandas as pd
from exhaustive.exhaustive.utils.utils import get_xtals_from_db

# from giant.jiffies.score_model import run as score_model
# from giant.jiffies.score_model import master_phil as score_phil
# import libtbx.phil

# Changing project just requires changing the directory.
# Note the /mnt prefix for running locally on laptop

#ini_folder = "/mnt/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios"
#ini_folder = "/mnt/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/titration_series"

# This works in ccp4 python. Currently using the score model script as the
# score params stuff is only needed if changing to a different residue
# such as GDP for janine

params.input.database_path = "/dls/labxchem/data/2016/lb13385-64/processing/database/soakDBDataFile.sqlite""

for xtal_name, pdb, mtz in get_xtals_from_db(params,
                                             refinement_outcomes="'4 - CompChem ready', "
                                                                 "'5 - Deposition ready',"
                                                                 "'6 - Deposited'" ):

    os.system("giant.score_model {} {}".format(pdb, mtz)))

# Currently pandas is failing to import in ccp4 python so this is done seperately in a conda env

# dfs = []
#
# for folder in os.listdir(ini_folder):
#     if os.path.exists(os.path.join(ini_folder,folder,"refine.pdb")):
#         edstats_csv = os.path.join(ini_folder, folder, "residue_scores.csv")
#         edstats_df = pd.read_csv(edstats_csv)
#         edstats_df['Dataset'] = folder
#         dfs.append(edstats_df)
#         print(folder)
#         print(edstats_df)
#         print(dfs)
#
#
# print("-------------------------------------------")
# compound_edstats = pd.concat(dfs, ignore_index=True)
# print(compound_edstats)
# compound_edstats.to_csv(os.path.join(ini_folder,'edstats.csv'))