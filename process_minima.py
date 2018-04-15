import os

import numpy as np

from process_exhaustive_search import write_minima_pdb
from validation import refine_after_exhasutive_search

# Initial Test

# working_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/NUDT7_Copied_atoms/NUDT7A-x1787/"
# input_pdb = os.path.join(working_dir,"multi-state-model.pdb")
# csv_name = os.path.join(working_dir,"u_iso_occupancy_vary_new_atoms")
# output_pdb = os.path.join(working_dir,"exhaustive_search_minima.pdb")
# #write_minima_pdb(input_pdb, output_pdb, csv_name)
#
# input_mtz = os.path.join(working_dir,"NUDT7A-x1787.free.mtz")
# input_cif = os.path.join(working_dir,"OX-210.cif")
# refine_params = os.path.join(working_dir,"multi-state-restraints.refmac.params")
# minima_pdb = os.path.join(working_dir,"exhaustive_search_minima.pdb")
#
# refine_after_exhasutive_search(input_pdb=minima_pdb, input_mtz=input_mtz, input_cif=input_cif,
#                                refine_params=refine_params, dataset_prefix="NUDT7A-x1787", working_dir=working_dir)

# Test with simulated_data

validation_path ="/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/validation_bound_ground/"

dataset_prefix = "NUDT7A-x1740"
folder_prefix = "NUDT7A-x1740_refine_occ_"
set_b = 40

for simul_occ in np.arange(0.05,0.96,0.05):

    working_dir = os.path.join(validation_path, folder_prefix + str(simul_occ).replace(".","_"))
    csv_name = "occ_{}_b_{}_u_iso".format(str(simul_occ).replace(".","_"),str(set_b))
    input_pdb = os.path.join(working_dir,
                             "{}_refine_occ_{}_set_b_{}.pdb".format(dataset_prefix,
                                                                    str(simul_occ).replace(".","_")),
                                                                    str(set_b))
    minima_pdb = os.path.join(working_dir,"exhaustive_search_minima.pdb")
    write_minima_pdb(input_pdb=input_pdb, output_pdb=minima_pdb, csv_name=os.path.join(validation_path,csv_name))
    input_cif = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1740/NUOOA000181a.cif
    simul_mtz = os.path.join(working_dir,"{}_simul_{}.mtz".format(dataset_prefix,
                                                                    str(simul_occ).replace(".","_")))
    refine_params = os.path.join(working_dir, "multi-state-restraints.refmac.params")

    refine_after_exhasutive_search(input_pdb=minima_pdb, input_mtz=simul_mtz, input_cif=input_cif,
                                   refine_params=refine_params, dataset_prefix=dataset_prefix, working_dir=working_dir)
#
# NUDT7A-x1740_refine_occ_0_05_set_b_40.pdb
# NUDT7A-x1740_simul_0_05.mtz