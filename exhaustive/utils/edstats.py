import os
import pandas as pd
# from giant.jiffies.score_model import run as score_model
# from giant.jiffies.score_model import master_phil as score_phil
# import libtbx.phil

# Changing project just requires changing the directory.
# Note the /mnt prefix for running locally on laptop

#ini_folder = "/mnt/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios"
ini_folder = "/mnt/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/titration_series"

# This works in ccp4 python. Currently using the score model script as the
# score params stuff is only needed if changing to a different residue
# such as GDP for janine

# for folder in os.listdir(ini_folder):
#     if os.path.exists(os.path.join(ini_folder,folder,"refine.pdb")):
# #         print(os.path.join(ini_folder,folder,"refine.pdb"))
# #         #score_params = score_phil.extract()
# #         #score_params.input.pdb1 = os.path.join(ini_folder,folder,"refine.pdb")
# #         #score_params.input.mtz1 = os.path.join(ini_folder,folder,"refine.mtz")
# #         #score_params.output.out_dir = os.path.join(ini_folder,folder,"edstats_on_refine")
# #         #score_params.selection.res_names= 'LIG'
# #         #score_model(score_params)
#         os.chdir(os.path.join(ini_folder,folder))
#         print(os.getcwd())
#         os.system("giant.score_model {} {}".format(os.path.join(ini_folder,folder,"refine.pdb"),
#                                                    os.path.join(ini_folder,folder,"refine.mtz")))

# Currently pandas is failing to import in ccp4 python so this is done seperately in a conda env

dfs = []

for folder in os.listdir(ini_folder):
    if os.path.exists(os.path.join(ini_folder,folder,"refine.pdb")):
        edstats_csv = os.path.join(ini_folder, folder, "residue_scores.csv")
        edstats_df = pd.read_csv(edstats_csv)
        edstats_df['Dataset'] = folder
        dfs.append(edstats_df)
        print(folder)
        print(edstats_df)
        print(dfs)


print("-------------------------------------------")
compound_edstats = pd.concat(dfs, ignore_index=True)
print(compound_edstats)
compound_edstats.to_csv(os.path.join(ini_folder,'edstats.csv'))