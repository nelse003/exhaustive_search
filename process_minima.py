import os

import numpy as np

# from process_exhaustive_search import write_minima_pdb
from process_exhaustive_search import get_minimum_fofc

#from validation import refine_after_exhasutive_search

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


set_b = 40
out_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/validation_bound_ground"

for simul_occ in np.arange(0.05,0.96,0.05):

    csv_name = "occ_{}_b_{}_u_iso".format(str(simul_occ).replace(".","_"),str(set_b))
    occ, u_iso, fo_fc = get_minimum_fofc(os.path.join(out_path,csv_name), b_fac = set_b)
    print("simul_occ {}, Occ {}, u_iso {} mean(|Fo -Fc| {}".format(simul_occ, occ,u_iso,fo_fc))

# Test with simulated_data

# validation_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/exhaustive_search_of_refmac_0_cyc/"
#
# dataset_prefix = "NUDT7A-x1740"
# folder_prefix = "NUDT7A-x1740_refine_occ_"
# set_b = 40
#
# for simul_occ in np.arange(0.05,0.96,0.05):
#
#     working_dir = os.path.join(validation_path, folder_prefix + str(simul_occ).replace(".","_"))
#     csv_name = "occ_{}_b_{}_u_iso".format(str(simul_occ).replace(".","_"),str(set_b))
#     input_pdb = os.path.join(working_dir,
#                              "{}_refine_occ_{}_set_b_{}.pdb".format(dataset_prefix,
#                                                                     str(simul_occ).replace(".","_"),
#                                                                     str(set_b)))
#     minima_pdb = os.path.join(working_dir,"exhaustive_search_minima.pdb")
#     write_minima_pdb(input_pdb=input_pdb, output_pdb=minima_pdb, csv_name=os.path.join(validation_path,csv_name))
#     input_cif = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1740/NUOOA000181a.cif"
#     simul_mtz = os.path.join(working_dir,"{}_simul_{}.mtz".format(dataset_prefix,
#                                                                     str(simul_occ).replace(".","_")))
#     refine_params = os.path.join(working_dir, "multi-state-restraints.refmac.params")
#
#     refine_after_exhasutive_search(input_pdb=minima_pdb, input_mtz=simul_mtz, input_cif=input_cif,
#                                    refine_params=refine_params, dataset_prefix=dataset_prefix, working_dir=working_dir)
#########################

# plot_random_refinement_starts(start_occ=0.05, end_occ=0.95, step=0.05,
#                               dataset_prefix=dataset_prefix, set_b=40, out_path=out_path)