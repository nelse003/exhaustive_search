import os
import random
import iotbx.mtz
import numpy as np
import libtbx.phil

# Local imports

from ..utils.refinement import refmac_0_cyc
from ..utils.utils import set_b_fac_all_occupancy_groups, wait_for_file_existence, get_csv_filepath, \
    set_b_fac_all_atoms, get_random_starting_occ_from_folder_name
from ..plotting.plot import scatter_plot, plot_3d_fofc_occ
from ..phil import master_phil, prepare_validate_phil, check_input_files

##############################################################
# Logging

##############################################################

def check_validate_input_files(params):

    """Check existence of input files, including all needed for validation code"""

    check_input_files(params)

    assert os.path.exists(params.output.out_dir), "{} does not exist".format(params.output.params.output.out_dir)
    assert os.path.exists(params.validate.input.params.validate.ground_state_pdb_path), \
        "Ground state pdb: \n{}\n does not exist".format(params.validate.input.params.validate.ground_state_pdb_path)
    assert os.path.exists(params.validate.input.params.validate.bound_state_pdb_path),\
        "Bound state pdb: \n{}\n does not exist".format(params.validate.input.params.validate.ground_state_pdb_path)

def occ_loop_merge_confs_simulate(params):

    """ Simulate Experimental data using phenix f_model. Run exhaustive_search on simulated data. 
    
    Loop over all occupancies between params.validate.options.start_simul_occ
    and params.validate.options.end_simul_occ, in sizes of params.validate.options.step_simulation.
     For each of these occupancies:
     > Run giant.merge_conformations to generate suitable (with occupancies set to (1- lig_occ) 
     for ground and to (lig_occ) for bound state) pdb file (multi-state-model.pdb) 
     to be passed to the simulation routine.
     > If a B factor is to be set (using params.validate.options.set_b) then set B factor of ground and bound states to 
     a fixed value using set_b_fac_all_occupancy_groups()
     > Simulate Fobs data using phenix.f_model, base ouput on reflections on params.input.mtz
     > Run exhaustive search routine on simulated data. Via qsub submission
     > Run phenix maps to get viewable map from simluated mtz.
    """

    check_validate_input_files(params = params)

    for lig_occupancy in np.arange(params.validate.options.start_simul_occ, 
                                   params.validate.options.end_simul_occ+params.validate.options.step_simulation/5, 
                                   params.validate.options.step_simulation):

        merged_pdb = os.path.join(params.output.out_dir,
                                  "{}_refine_occ_{}.pdb".format(params.input.xtal_name, 
                                                                str(lig_occupancy).replace(".", "_")))

        if params.exhaustive.options.overwrite or not os.path.exists(os.path.join(merged_pdb)):
            os.system("giant.merge_conformations input.major={} input.minor={} "
                      "major_occupancy={} minor_occupancy={} output.pdb={}".format(
                params.validate.ground_state_pdb_path, params.validate.bound_state_pdb_path,
                str(1 - lig_occupancy), str(lig_occupancy), merged_pdb))

        if params.validate.options.set_b is not None:
            merged_file_name, _ = os.path.splitext(merged_pdb)

            if params.exhaustive.options.overwrite or not os.path.exists(merged_pdb
             + "_set_b_{}.pdb".format(str(params.validate.options.set_b).replace(".", "_"))):

                # set_b_fac_all_atoms(input_pdb = merged_pdb,
                #                     output_pdb = merged_file_name + "_set_all_b_{}.pdb".format(
                #                        str(params.validate.options.set_b).replace(".", "_")),
                #                     b_fac = params.validate.options.set_b)

                set_b_fac_all_occupancy_groups(input_pdb = merged_pdb,
                                               output_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                                   str(params.validate.options.set_b).replace(".", "_")),
                                               b_fac = params.validate.options.set_b)

            merged_pdb = merged_file_name + "_set_b_{}.pdb".format(str(params.validate.options.set_b).replace(".", "_"))

        os.chdir(params.output.out_dir)

        if params.exhaustive.options.overwrite or not \
                os.path.exists(os.path.join(params.output.out_dir, merged_pdb +".mtz")):

            #TODO Work out data column label for sensible input?

            if merged_pdb != merged_file_name + "_set_b_{}.pdb".format(str(params.validate.options.set_b).replace(".","_")):
                exit()

            os.system("phenix.fmodel data_column_label=\"F,SIGF\" {} {} type=real".format(merged_pdb, params.input.mtz))
            #os.system("phenix.fmodel high_res={} type=real {}".format(high, merged_pdb, merged_pdb +".mtz" ))

        assert os.path.exists(merged_pdb+".mtz")

        sh_file = "{}_occ_{}_b_{}.sh".format(params.input.xtal_name,
                                             str(lig_occupancy).replace(".", "_"),
                                             str(params.validate.options.set_b).replace(".", "_"))


        if params.exhaustive.options.overwrite or not os.path.exists(os.path.join(params.output.out_dir,sh_file)):
            with open(os.path.join(params.output.out_dir, sh_file),'w') as file:

                file.write("#!/bin/bash\n")
                file.write("export XChemExplorer_DIR=\"/dls/science/"
                           "groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
                file.write("source /dls/science/groups/i04-1/software/"
                           "XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")

                file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/"
                           "elliot-dev/Work/exhaustive_search/exhaustive/exhaustive.py"
                           "input.pdb={} input.mtz={} output.out_dir={} xtal_name={} "
                           "options.csv_name={} options.step={} options.buffer={} "
                           "options.params.exhaustive.options.grid_spacing={} "
                           "params.exhaustive.options.generate_mtz={}".format(merged_pdb, merged_pdb +".mtz",
                                                                            params.output.out_dir, 
                                                                            params.input.xtal_name, 
                                                                            params.exhaustive.options.csv_name,
                                                                            params.exhaustive.options.step,
                                                                            params.validate.options.buffer,
                                                                            params.exhaustive.options.grid_spacing,
                                                                            params.exhaustive.options.generate_mtz))


        if params.exhaustive.options.overwrite or not os.path.exists(os.path.join(
                params.output.out_dir,params.exhaustive.csv_name +".csv")):
            os.system("qsub -o {} -e {} {}".format(os.path.join(params.output.out_dir,"output_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                                   os.path.join(params.output.out_dir,"error_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                                   os.path.join(params.output.out_dir, sh_file)))

        if params.exhaustive.options.overwrite or not os.path.exists(merged_pdb +".mtz") or not os.path.exists(merged_pdb):
            os.system("phenix.maps {} {} maps.map.map_type=\"mfo-Dfc\"".format(merged_pdb, merged_pdb +".mtz"))

def run(params):

    params = master_phil.extract()
    # params.input.in_path = "/dls/labxchem/data/2016/lb13385-61/processing/analysis/initial_model/FALZA-x0085"
    # params.input.mtz = os.path.join(params.input.in_path, "FALZA-x0085.free.mtz")
    modified_phil = master_phil.format(python_object=params)
    modified_phil = prepare_validate_phil(modified_phil)
    params = modified_phil.extract()
    check_validate_input_files(params)

    # in_path = "/dls/labxchem/data/2016/lb13385-61/processing/analysis/initial_model/FALZA-x0085"
    # params.validate.bound_state_pdb_path = os.path.join(in_path,"refine.output.bound-state.pdb")
    # params.validate.ground_state_pdb_path =  os.path.join(in_path,"refine.output.ground-state.pdb")
    # input_cif = os.path.join(in_path, "FMOPL000287a.cif")
    # params.input.xtal_name = "FALZA-x0085"
    # params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation.py/exhaustive_search_phenix_fmodel/FALZA-x0085-code-cleanup-tests"
    # params.validate.options.set_b= 40

    if not os.path.exists(params.output.out_dir):
        os.mkdir(params.output.out_dir)

    # This loop runs exhaustive search many times across simulated data
    occ_loop_merge_confs_simulate(params)

    # Waits for occupancy csvs to be output
    for file_path in get_csv_filepath(params.output.out_dir,
                                      set_b=params.validate.options.set_b,
                                      step= params.validate.,
                                      start_occ=0.05,
                                      end_occ=0.95):
        wait_for_file_existence(file_path, wait_time=10000)

    # This plots exhaustive search results, to confirm whether exhaustive search recovers the simulated occupancy
    os.chdir(params.output.out_dir)
    plot_3d_fofc_occ(0.05, 0.95, step=0.05, set_b=40, xtal_name=params.input.xtal_name)


    os.chdir(params.output.out_dir)
    for simul_occ in np.arange(0.05, 0.95, 0.05):
        csv_name = "occ_{}_b_{}_u_iso".format(str(simul_occ).replace(".", "_"),params.validate.options.set_b)
        scatter_plot(csv_name, title_text="Phenix.fmodel at occ {}".format(simul_occ))

if(__name__ == "__main__"):
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args = sys.argv[1:])

# ################################################
# # Check and turn to functions
# ################################################
#
#
# # occ_loop_merge_confs_simulate_with_refmac_0(params.validate.bound_state_pdb_path,
# #                                          params.validate.ground_state_pdb_path,
# #                                          params.input.mtz,
# #                                          params.input.xtal_name,
# #                                          params.output.out_dir,
# #                                          params.validate.options.set_b = 40)
# # submit_exhasutive_with_refmac_0(params.input.xtal_name, params.output.out_dir, params.validate.options.set_b = 40)
# # os.chdir(params.output.out_dir)
# # plot_fofc_occ(0.05, 0.95, 0.05, 40)
#
# validation_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation.py/validation_bound_ground/"
#
# params.input.xtal_name = "NUDT7A-x1740"
# folder_prefix = "NUDT7A-x1740_refine_occ_"
# params.validate.options.set_b = 40
#
# for simul_occ in np.arange(0.05, 0.96, 0.05):
#
#     working_dir = os.path.join(validation_path, folder_prefix + str(simul_occ).replace(".", "_"))
#
#     # Create a folder for each
#     # giant.score_model for simulated data: use refmac 0 cycles version for compatibility
#
#     # if not os.path.exists(os.path.join(working_dir, "simulated_refmac_0_score_model")):
#     #     os.mkdir(os.path.join(working_dir, "simulated_refmac_0_score_model"))
#     # os.chdir(os.path.join(working_dir, "simulated_refmac_0_score_model"))
#     #
#     # input_pdb = os.path.join(working_dir,
#     #                          "{}_occ_{}_b_{}_refmac_0_cyc.pdb".format(params.input.xtal_name, str(simul_occ).replace(".","_"),
#     #                                                                             str(params.validate.options.set_b).replace(".","_")))
#     # params.input.mtz = os.path.join(working_dir,
#     #                          "{}_occ_{}_b_{}_refmac_0_cyc.mtz".format(params.input.xtal_name, str(simul_occ).replace(".","_"),
#     #                                                                             str(params.validate.options.set_b).replace(".","_")))
#     # os.system("giant.score_model input.pdb1={} input.mtz1={}".format(input_pdb,params.input.mtz))
#     #
#     # # giant.score_model for exhaustive search minima pdb
#     #
#     # if not os.path.exists(os.path.join(working_dir, "exhaustive_search_minima_score_model")):
#     #     os.mkdir(os.path.join(working_dir, "exhaustive_search_minima_score_model"))
#     # os.chdir(os.path.join(working_dir, "exhaustive_search_minima_score_model"))
#     #
#     # input_pdb = os.path.join(working_dir,
#     #                          "exhaustive_seach_minima.pdb".format(str(simul_occ).replace(".","_")))
#     #
#     # os.system("giant.score_model input.pdb1={} input.mtz1={}".format(input_pdb,params.input.mtz))
#     #
#     # # giant.score_model for exhaustive search minima, after refinement
#     #
#     # if not os.path.exists(os.path.join(working_dir, "exhaustive_search_minima_refined_score_model")):
#     #     os.mkdir(os.path.join(working_dir, "exhaustive_search_minima_refined_score_model"))
#     # os.chdir(os.path.join(working_dir, "exhaustive_search_minima_refined_score_model"))
#     #
#     # folders = [name for name in os.listdir(working_dir) if os.path.isdir(os.path.join(working_dir, name))]
#     #
#     # es_refine = []
#     # for folder in folders:
#     #     if folder.find('refine_after_exhaustive_search') != -1:
#     #         es_refine.append(int(folder[-4:]))
#     # es_refine_folder = os.path.join(ES_folder,
#     #                                 "{}_refine_after_exhaustive_search{}".format(params.input.xtal_name,
#     #                                                                              str(max(es_refine)).rjust(4, '0')))
#     # ES_refine_pdb = os.path.join(es_refine_folder, "{}_refine_after_exhaustive_search.pdb".format(params.input.xtal_name))
#     # ES_refine_mtz = os.path.join(es_refine_folder, "{}_refine_after_exhaustive_search.mtz".format(params.input.xtal_name))
#     #
#     # os.system("giant.score_model input.pdb1={} input.mtz1={} input.pdb2={} input.mtz2={}".format(ES_refine_pdb,
#     #                                                                                              ES_refine_mtz,
#     #                                                                                              input_pdb,
#     #                                                                                              params.input.mtz))
#
#     # giant.score model for refined from random point (5 cycles)
#
#     for starting_rand_occ in get_random_starting_occ_from_folder_name(simul_occ, validation_path, params.input.xtal_name):
#         cur_dir = os.path.join(working_dir,
#                                params.input.xtal_name + "_expected_occ_"
#                                + str(simul_occ).replace(".", "_") + "_b_" + str(params.validate.options.set_b) + "_supplied_occ_" +
#                                str(starting_rand_occ).replace(".", "_"))
#
#         print(cur_dir)
#
#     # giant.score_model for refined from random point (20 cycles)
#
#
#     exit()
#
#
#     # plot_random_refinement_with_ES(start_occ=0.05, end_occ=0.95, step=0.05,
#     #                               params.input.xtal_name=params.input.xtal_name, params.validate.options.set_b=40, params.output.out_dir=params.output.out_dir)
#
#     # quick_refine_repeats(start_occ = 0.05, end_occ = 0.95, step = 0.05,
#     #                      params.input.xtal_name = params.input.xtal_name, params.validate.options.set_b=40, params.output.out_dir = params.output.out_dir, input_cif = input_cif)
#
#
#     # occ_loop_merge_refine_random_confs_simulate(params.validate.bound_state_pdb_path,
#     #                                              params.validate.ground_state_pdb_path,
#     #                                              params.input.mtz,
#     #                                              params.input.xtal_name,
#     #                                              params.output.out_dir,
#     #                                              input_cif,
#     #                                              params.validate.options.set_b = 40)
