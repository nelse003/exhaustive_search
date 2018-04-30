import os
import random

import iotbx.mtz
import numpy as np

from plot_exhaustive_search import scatter_plot, plot_fofc_occ
from refinement import refmac_0_cyc
from utils import set_u_iso_all_occupancy_groups, wait_for_file_existence, get_csv_filepath


def occ_loop_merge_confs_simulate(bound_state_pdb_path,
                                  ground_state_pdb_path,
                                  input_mtz,
                                  dataset_prefix,
                                  out_path,
                                  set_b = None,
                                  step_simul = 0.05,
                                  step_sampling = 0.01,
                                  start_simul_occ = 0.05,
                                  end_simul_occ = 0.95,
                                  buffer = 0,
                                  grid_spacing = 0.25,
                                  overwrite = False):



    for lig_occupancy in np.arange(start_simul_occ, end_simul_occ+step_simul/5, step_simul):

        csv_name = "occ_{}_b_{}_u_iso".format(str(lig_occupancy).replace(".", "_"), str(set_b).replace(".", "_"))
        # Skip if already generated csv
        if not overwrite:
            if os.path.exists(os.path.join(out_path,csv_name +".csv")):
                print("Skipping {} because {} exists".format(lig_occupancy,csv_name) )
                continue

        merged_pdb = os.path.join(out_path,
                                  "{}_refine_occ_{}.pdb".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        os.system("giant.merge_conformations input.major={} input.minor={} "
                  "major_occupancy={} minor_occupancy={} output.pdb={}".format(
            ground_state_pdb_path, bound_state_pdb_path, str(1 - lig_occupancy), str(lig_occupancy), merged_pdb))

        simulate_log = os.path.join(out_path,"{}_simul_{}.log".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))
        simulate_mtz = os.path.join(out_path,"{}_simul_{}.mtz".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        if set_b is not None:
            merged_file_name, _ = os.path.splitext(merged_pdb)

            set_u_iso_all_occupancy_groups(input_pdb = merged_pdb,
                                           output_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_")),
                                           b_fac = set_b)

            merged_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_"))

        # os.system("ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/simulate_experimental_data.py "
        #           "input.xray_data.file_name={} "
        #           "model.file_name={} input.xray_data.label=\"F,SIGF\" "
        #           "output.logfile={} output.hklout={}".format(input_mtz, merged_pdb,
        #                                                       simulate_log, simulate_mtz))

        os.chdir(out_path)

        o = iotbx.mtz.object(input_mtz)
        low,high =o.max_min_resolution()

        os.system("phenix.fmodel high_res={} type=real {} ".format(high, merged_pdb))

        # Exhaustive search
        sh_file = "{}_occ_{}_b_{}.sh".format(dataset_prefix,
                                             str(lig_occupancy).replace(".", "_"),
                                             str(set_b).replace(".", "_"))

        with open(os.path.join(out_path, sh_file),'w') as file:

            file.write("#!/bin/bash\n")
            file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
            file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")

            file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive_search.py"
                       " input.pdb={} input.mtz={} output.out_dir={} xtal_name={} "
                       "options.csv_name={} options.step={} options.buffer={} "
                       "options.grid_spacing={}".format(merged_pdb, merged_pdb +".mtz", out_path, dataset_prefix, csv_name,
                                                        step_sampling, buffer, grid_spacing))

        # os.chmod(os.path.join(out_path, sh_file),0777)
        # os.system(os.path.join(out_path, sh_file))
        # exit()

        os.system("qsub -o {} -e {} {}".format(os.path.join(out_path,"output_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                               os.path.join(out_path,"error_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                               os.path.join(out_path, sh_file)))

        #Refmac 0 cycles
        output_pdb = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.pdb".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        output_mtz = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.mtz".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        refmac_0_cyc(input_mtz = simulate_mtz, input_pdb = merged_pdb,
                     output_pdb = output_pdb , output_mtz = output_mtz,
                     occupancy= lig_occupancy)

def occ_loop_merge_confs_simulate_with_refmac_0(bound_state_pdb_path,
                                             ground_state_pdb_path,
                                             input_mtz,
                                             dataset_prefix,
                                             out_path,
                                             set_b = None):

    """ Run exhaustive search on files after refmac 0 cycles"""

    for lig_occupancy in np.arange(0.05, 0.96, 0.05):
        merged_pdb = os.path.join(out_path,
                                  "{}_refine_occ_{}.pdb".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        os.system("giant.merge_conformations input.major={} input.minor={} "
                  "major_occupancy={} minor_occupancy={} output.pdb={}".format(
            ground_state_pdb_path, bound_state_pdb_path, str(1 - lig_occupancy), str(lig_occupancy), merged_pdb))

        simulate_log = os.path.join(out_path,"{}_simul_{}.log".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))
        simulate_mtz = os.path.join(out_path,"{}_simul_{}.mtz".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        if set_b is not None:
            merged_file_name, _ = os.path.splitext(merged_pdb)

            set_u_iso_all_occupancy_groups(input_pdb = merged_pdb,
                                           output_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_")),
                                           b_fac = set_b)

            merged_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_"))

        os.system("ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/simulate_experimental_data.py "
                  "input.xray_data.file_name={} "
                  "model.file_name={} input.xray_data.label=\"F,SIGF\" "
                  "output.logfile={} output.hklout={}".format(input_mtz, merged_pdb, simulate_log, simulate_mtz))

        # Exhaustive search
        sh_file = "{}_occ_{}_b_{}.sh".format(dataset_prefix,
                                             str(lig_occupancy).replace(".", "_"),
                                             str(set_b).replace(".", "_"))
        csv_name = "occ_{}_b_{}_u_iso".format(str(lig_occupancy).replace(".", "_"),str(set_b).replace(".", "_"))

        #Refmac 0 cycles
        output_pdb = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.pdb".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        output_mtz = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.mtz".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        refmac_0_cyc(input_mtz = simulate_mtz, input_pdb = merged_pdb,
                     output_pdb = output_pdb , output_mtz = output_mtz,
                     occupancy= lig_occupancy)


        with open(os.path.join(out_path, sh_file),'w') as file:

            file.write("#!/bin/bash\n")
            file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
            file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")

            file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive_search.py"
                       " input.pdb={} input.mtz={} output.out_dir={} xtal_name={} "
                       "options.csv_name={}".format(output_pdb,output_mtz,out_path,dataset_prefix,csv_name))


def occ_loop_merge_refine_random_confs_simulate(bound_state_pdb_path,
                                             ground_state_pdb_path,
                                             input_mtz,
                                             dataset_prefix,
                                             out_path,
                                             input_cif,
                                             set_b = None):

    for lig_occupancy in np.arange(0.05, 0.96, 0.05):


        out_path = os.path.join(out_path, "{}_refine_occ_{}".format(dataset_prefix,
                                                                    str(lig_occupancy).replace(".", "_")))
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        os.chdir(out_path)

        merged_pdb = os.path.join(out_path,
                                  "{}_refine_occ_{}.pdb".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        os.system("giant.merge_conformations input.major={} input.minor={} "
                  "major_occupancy={} minor_occupancy={} output.pdb={}".format(
            ground_state_pdb_path, bound_state_pdb_path, str(1 - lig_occupancy), str(lig_occupancy), merged_pdb))

        simulate_log = os.path.join(out_path,"{}_simul_{}.log".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))
        simulate_mtz = os.path.join(out_path,"{}_simul_{}.mtz".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        if set_b is not None:
            merged_file_name, _ = os.path.splitext(merged_pdb)

            set_u_iso_all_occupancy_groups(input_pdb = merged_pdb,
                                           output_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_")),
                                           b_fac = set_b)

            merged_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_"))

        os.system("ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/simulate_experimental_data.py "
                  "input.xray_data.file_name={} "
                  "model.file_name={} input.xray_data.label=\"F,SIGF\" "
                  "output.logfile={} output.hklout={}".format(input_mtz, merged_pdb, simulate_log, simulate_mtz))


        num_random_starts = 0
        while num_random_starts < 30:

            num_random_starts += 1
            # Make a pdb file with the occupancy of the ligand set to the random value

            refinement_starting_occ = random.random()

            out_path = os.path.join(out_path,"{}_expected_occ_{}_b_{}_supplied_occ_{}".format(dataset_prefix,
                                                                      str(lig_occupancy).replace(".", "_"),
                                                                      str(set_b).replace(".", "_"),
                                                                      str(refinement_starting_occ).replace(".", "_")))
            if not os.path.exists(out_path):
                os.mkdir(out_path)
            os.chdir(out_path)

            refinement_random_occ_pdb = os.path.join(out_path,"{}_random_refine_occ_{}.pdb".format(
                                                                    dataset_prefix,
                                                                    str(refinement_starting_occ).replace(".", "_")))

            os.system("giant.merge_conformations input.major={} input.minor={} "
                      "major_occupancy={} minor_occupancy={} output.pdb={}".format(
                        ground_state_pdb_path,
                        bound_state_pdb_path,
                        str(1 - refinement_starting_occ),
                        str(refinement_starting_occ),
                        refinement_random_occ_pdb))

            out_prefix = "expected_occ_{}_supplied_occ_{}".format(str(lig_occupancy).replace(".", "_"),
                                                                 str(refinement_starting_occ).replace(".", "_"))
            dir_prefix = "refine_" + out_prefix +"_"
            params = "multi-state-restraints.refmac.params"

            sh_file = "{}_expected_occ_{}_b_{}_supplied_occ_{}.sh".format(dataset_prefix,
                                                 str(lig_occupancy).replace(".", "_"),
                                                 str(set_b).replace(".", "_"),
                                                 str(refinement_starting_occ).replace(".", "_"))

            with open(os.path.join(out_path, sh_file),'w') as file:

                file.write("#!/bin/bash\n")
                file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
                file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
                file.write("giant.quick_refine input.pdb={} input.mtz={} input.cif={} "
                           "input.params={} output.out_prefix={} output.dir_prefix={}".format(
                            refinement_random_occ_pdb, simulate_mtz, input_cif,params, out_prefix, dir_prefix))

            os.system("qsub -o {} -e {} {}".format(os.path.join(out_path,"output_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                                   os.path.join(out_path,"error_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                                   os.path.join(out_path, sh_file)))

            out_path = os.path.dirname(out_path)

        #Refmac 0 cycles
        output_pdb = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.pdb".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        output_mtz = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.mtz".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        refmac_0_cyc(input_mtz = simulate_mtz, input_pdb = merged_pdb,
                     output_pdb = output_pdb , output_mtz = output_mtz,
                     occupancy= lig_occupancy)

        out_path = os.path.dirname(out_path)
        print(out_path)


in_path = "/dls/labxchem/data/2018/lb18145-55/processing/analysis/initial_model/NUDT22A-x1058"
bound_state_pdb_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/NUDT22A_new_ground_state_x1058/multi-state-model.split.bound-state.pdb"
ground_state_pdb_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/NUDT22A_new_ground_state_x1058/multi-state-model.split.ground-state.pdb"
input_mtz = os.path.join(in_path, "NUDT22A-x1058.free.mtz")
# input_cif = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1740/NUOOA000181a.cif"
dataset_prefix = "NUDT22A-x1058"
out_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/exhaustive_search_phenix_fmodel/NUDT22A-x1058"

if not os.path.exists(out_path):
    os.mkdir(out_path)

# This loop runs exhaustive search many times across simulated data
occ_loop_merge_confs_simulate(bound_state_pdb_path,
                              ground_state_pdb_path,
                              input_mtz,
                              dataset_prefix,
                              out_path,
                              set_b = 40,
                              step_simul= 0.01,
                              start_simul_occ= 0.01,
                              end_simul_occ= 0.99,
                              buffer = 0,
                              grid_spacing = 0.25)


# Waits for occupancy csvs to be output
for file_path in get_csv_filepath(out_path, set_b=40, step=0.01, start_occ=0.01, end_occ=0.99):
    wait_for_file_existence(file_path, wait_time=1000)

# This plots exhaustive search results, to confirm whether exhaustive search recovers the simulated occupancy
os.chdir(out_path)
plot_fofc_occ(0.01, 0.99, step=0.01, set_b=40, dataset_prefix=dataset_prefix)


os.chdir(out_path)
for simul_occ in np.arange(0.01, 0.99, 0.01):
    csv_name = "occ_{}_b_40_u_iso".format(str(simul_occ).replace(".", "_"))
    scatter_plot(csv_name, title_text="Phenix.fmodel at occ {}".format(simul_occ))
exit()

################################################
# Check and turn to functions
################################################


# occ_loop_merge_confs_simulate_with_refmac_0(bound_state_pdb_path,
#                                          ground_state_pdb_path,
#                                          input_mtz,
#                                          dataset_prefix,
#                                          out_path,
#                                          set_b = 40)
# submit_exhasutive_with_refmac_0(dataset_prefix, out_path, set_b = 40)
# os.chdir(out_path)
# plot_fofc_occ(0.05, 0.95, 0.05, 40)



validation_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/validation_bound_ground/"

dataset_prefix = "NUDT7A-x1740"
folder_prefix = "NUDT7A-x1740_refine_occ_"
set_b = 40

for simul_occ in np.arange(0.05, 0.96, 0.05):

    working_dir = os.path.join(validation_path, folder_prefix + str(simul_occ).replace(".", "_"))

    # Create a folder for each
    # giant.score_model for simulated data: use refmac 0 cycles version for compatibility

    # if not os.path.exists(os.path.join(working_dir, "simulated_refmac_0_score_model")):
    #     os.mkdir(os.path.join(working_dir, "simulated_refmac_0_score_model"))
    # os.chdir(os.path.join(working_dir, "simulated_refmac_0_score_model"))
    #
    # input_pdb = os.path.join(working_dir,
    #                          "{}_occ_{}_b_{}_refmac_0_cyc.pdb".format(dataset_prefix, str(simul_occ).replace(".","_"),
    #                                                                             str(set_b).replace(".","_")))
    # input_mtz = os.path.join(working_dir,
    #                          "{}_occ_{}_b_{}_refmac_0_cyc.mtz".format(dataset_prefix, str(simul_occ).replace(".","_"),
    #                                                                             str(set_b).replace(".","_")))
    # os.system("giant.score_model input.pdb1={} input.mtz1={}".format(input_pdb,input_mtz))
    #
    # # giant.score_model for exhaustive search minima pdb
    #
    # if not os.path.exists(os.path.join(working_dir, "exhaustive_search_minima_score_model")):
    #     os.mkdir(os.path.join(working_dir, "exhaustive_search_minima_score_model"))
    # os.chdir(os.path.join(working_dir, "exhaustive_search_minima_score_model"))
    #
    # input_pdb = os.path.join(working_dir,
    #                          "exhaustive_seach_minima.pdb".format(str(simul_occ).replace(".","_")))
    #
    # os.system("giant.score_model input.pdb1={} input.mtz1={}".format(input_pdb,input_mtz))
    #
    # # giant.score_model for exhaustive search minima, after refinement
    #
    # if not os.path.exists(os.path.join(working_dir, "exhaustive_search_minima_refined_score_model")):
    #     os.mkdir(os.path.join(working_dir, "exhaustive_search_minima_refined_score_model"))
    # os.chdir(os.path.join(working_dir, "exhaustive_search_minima_refined_score_model"))
    #
    # folders = [name for name in os.listdir(working_dir) if os.path.isdir(os.path.join(working_dir, name))]
    #
    # es_refine = []
    # for folder in folders:
    #     if folder.find('refine_after_exhaustive_search') != -1:
    #         es_refine.append(int(folder[-4:]))
    # es_refine_folder = os.path.join(ES_folder,
    #                                 "{}_refine_after_exhaustive_search{}".format(dataset_prefix,
    #                                                                              str(max(es_refine)).rjust(4, '0')))
    # ES_refine_pdb = os.path.join(es_refine_folder, "{}_refine_after_exhaustive_search.pdb".format(dataset_prefix))
    # ES_refine_mtz = os.path.join(es_refine_folder, "{}_refine_after_exhaustive_search.mtz".format(dataset_prefix))
    #
    # os.system("giant.score_model input.pdb1={} input.mtz1={} input.pdb2={} input.mtz2={}".format(ES_refine_pdb,
    #                                                                                              ES_refine_mtz,
    #                                                                                              input_pdb,
    #                                                                                              input_mtz))

    # giant.score model for refined from random point (5 cycles)

    for starting_rand_occ in get_random_starting_occ_from_folder_name(simul_occ, validation_path, dataset_prefix):
        cur_dir = os.path.join(working_dir,
                               dataset_prefix + "_expected_occ_"
                               + str(simul_occ).replace(".", "_") + "_b_" + str(set_b) + "_supplied_occ_" +
                               str(starting_rand_occ).replace(".", "_"))

        print(cur_dir)

    # giant.score_model for refined from random point (20 cycles)


    exit()


    # plot_random_refinement_with_ES(start_occ=0.05, end_occ=0.95, step=0.05,
    #                               dataset_prefix=dataset_prefix, set_b=40, out_path=out_path)

    # quick_refine_repeats(start_occ = 0.05, end_occ = 0.95, step = 0.05,
    #                      dataset_prefix = dataset_prefix, set_b=40, out_path = out_path, input_cif = input_cif)


    # occ_loop_merge_refine_random_confs_simulate(bound_state_pdb_path,
    #                                              ground_state_pdb_path,
    #                                              input_mtz,
    #                                              dataset_prefix,
    #                                              out_path,
    #                                              input_cif,
    #                                              set_b = 40)
