import os
import random

import iotbx.mtz
import numpy as np

from plot_exhaustive_search import scatter_plot, plot_3d_fofc_occ
from refinement import refmac_0_cyc
from utils import set_b_fac_all_occupancy_groups, wait_for_file_existence, get_csv_filepath


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
                                  overwrite = False,
                                  input_cif = None):

    """ Simulate Experimental data using phenix f_model. Run exhaustive_search on simulated data. 
    
    Loop over all occupancies between start_simul_occ and end_simul_occ, in sizes of step_simul. 
     For each of these occupancies:
     > Run giant.merge_conformations to generate suitable (with occupancies set to (1- lig_occ) 
     for ground and to (lig_occ) for bound state) pdb file (multi-state-model.pdb) 
     to be passed to the simulation routine.
     > If a B factor is to be set (using set_b) then set B factor of ground and bound states to 
     a fixed value using set_u_iso_all_occupancy_groups()
     > Get high resolution shell from input_mtz, use the free.mtz. Needed for simulating data using 
     phenix.f_model
     > Simulate Fobs data using phenix.f_model
     > Run exhaustive search routine on simulated data. Via qsub submission
     > Run refmac with 0 cycles of refinement to get viewable mtz from siluated mtz.
    """

    assert os.path.exists(out_path), "{} does not exist".format(out_path)
    assert os.path.exists(ground_state_pdb_path), "Ground state pdb: \n{}\n does not exist".format(ground_state_pdb_path)
    assert os.path.exists(bound_state_pdb_path),"bound state pdb: \n{}\n does not exist".format(bound_state_pdb_path)
    assert os.path.exists(input_mtz), "Input mtz: \n{}\ndoes not exist".format(input_mtz)
    assert os.path.exists(input_cif), "Input cif:\n{}\ndoes not exist".format(input_cif)

    for lig_occupancy in np.arange(start_simul_occ, end_simul_occ+step_simul/5, step_simul):

        csv_name = "occ_{}_b_{}_u_iso".format(str(lig_occupancy).replace(".", "_"), str(set_b).replace(".", "_"))

        merged_pdb = os.path.join(out_path,
                                  "{}_refine_occ_{}.pdb".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        if overwrite or not os.path.exists(os.path.join(merged_pdb)):
            os.system("giant.merge_conformations input.major={} input.minor={} "
                      "major_occupancy={} minor_occupancy={} output.pdb={}".format(
                ground_state_pdb_path, bound_state_pdb_path, str(1 - lig_occupancy), str(lig_occupancy), merged_pdb))

        if set_b is not None:
            merged_file_name, _ = os.path.splitext(merged_pdb)

            set_b_fac_all_occupancy_groups(input_pdb = merged_pdb,
                                           output_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_")),
                                           b_fac = set_b)

            merged_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_"))

        simulate_mtz = os.path.join(out_path,
                                    "{}_simul_{}.mtz".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        simulate_log = os.path.join(out_path,"{}_simul_{}.log".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))
        # os.system("ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/simulate_experimental_data.py "
        #           "input.xray_data.file_name={} "
        #           "model.file_name={} input.xray_data.label=\"F,SIGF\" "
        #           "output.logfile={} output.hklout={}".format(input_mtz, merged_pdb,
        #                                                       simulate_log, simulate_mtz))

        os.chdir(out_path)

        if overwrite or not os.path.exists(os.path.join(out_path, merged_pdb +".mtz")):

            o = iotbx.mtz.object(input_mtz)
            low,high =o.max_min_resolution()
            #print("phenix.fmodel data_column_label=\"F,SIGF\" {} {} type=real".format(merged_pdb, input_mtz ))

            #TODO Work out data column label for sensible input?

            os.system("phenix.fmodel data_column_label=\"F,SIGF\" {} {} type=real".format(merged_pdb, input_mtz))
            #os.system("phenix.fmodel high_res={} type=real {}".format(high, merged_pdb, merged_pdb +".mtz" ))

        assert os.path.exists(merged_pdb+".mtz")

        sh_file = "{}_occ_{}_b_{}.sh".format(dataset_prefix,
                                             str(lig_occupancy).replace(".", "_"),
                                             str(set_b).replace(".", "_"))

        if lig_occupancy == 0.05:
            generate_mtz= True
        else:
            generate_mtz = False

        if overwrite or os.path.exists(os.path.join(out_path,sh_file)):
            with open(os.path.join(out_path, sh_file),'w') as file:

                file.write("#!/bin/bash\n")
                file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
                file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")

                file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive_search.py"
                           " input.pdb={} input.mtz={} output.out_dir={} xtal_name={} "
                           "options.csv_name={} options.step={} options.buffer={} "
                           "options.grid_spacing={} generate_mtz={}".format(merged_pdb, merged_pdb +".mtz", out_path, dataset_prefix, csv_name,
                                                            step_sampling, buffer, grid_spacing, generate_mtz))


        if overwrite or not os.path.exists(os.path.join(out_path,csv_name +".csv")):

            os.system("qsub -o {} -e {} {}".format(os.path.join(out_path,"output_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                                   os.path.join(out_path,"error_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                                   os.path.join(out_path, sh_file)))

        output_pdb = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.pdb".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        output_mtz = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.mtz".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        output_cif =os.path.join(out_path, os.path.basename(input_cif))

        if overwrite or not os.path.exists(output_mtz) or not os.path.exists(output_pdb):
            # refmac_0_cyc(input_mtz = simulate_mtz, input_pdb = merged_pdb,
            #              output_pdb = output_pdb , output_mtz = output_mtz,
            #              input_cif = input_cif,output_cif= output_cif,
            #              occupancy= lig_occupancy)
            # os.system("phenix.refine {} {} {} main.number_of_macro_cycles=0 "
            #           "refinement.input.xray_data.r_free_flags.generate=True".format(merged_pdb,
            #                                                                          merged_pdb +".mtz",
            #                                                                          input_cif))
            os.system("phenix.maps {} {} maps.map.map_type=\"mfo-Dfc\"".format(merged_pdb, merged_pdb +".mtz"))

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

            set_b_fac_all_occupancy_groups(input_pdb = merged_pdb,
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

            set_b_fac_all_occupancy_groups(input_pdb = merged_pdb,
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

def gradient(csv_name):
    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)
    occ = data[:,0]
    u_iso = data[:,2]
    fo_fc = data[:,3]



in_path = "/dls/labxchem/data/2016/lb13385-61/processing/analysis/initial_model/FALZA-x0085"
bound_state_pdb_path = os.path.join(in_path,"refine.output.bound-state.pdb")
ground_state_pdb_path =  os.path.join(in_path,"refine.output.ground-state.pdb")
input_mtz = os.path.join(in_path, "FALZA-x0085.free.mtz")
input_cif = os.path.join(in_path, "FMOPL000287a.cif")
dataset_prefix = "FALZA-x0085"
out_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/exhaustive_search_phenix_fmodel/FALZA-x0085"
set_b= 40

if not os.path.exists(out_path):
    os.mkdir(out_path)

# # This loop runs exhaustive search many times across simulated data
# occ_loop_merge_confs_simulate(bound_state_pdb_path,
#                               ground_state_pdb_path,
#                               input_mtz,
#                               dataset_prefix,
#                               out_path,
#                               set_b = set_b,
#                               step_simul= 0.05,
#                               start_simul_occ= 0.05,
#                               end_simul_occ= 0.95,
#                               buffer = 1,
#                               grid_spacing = 0.25,
#                               overwrite = False,
#                               input_cif = input_cif)
#
#
# # Waits for occupancy csvs to be output
# for file_path in get_csv_filepath(out_path, set_b=set_b, step=0.05, start_occ=0.05, end_occ=0.95):
#     wait_for_file_existence(file_path, wait_time=10000)

# This plots exhaustive search results, to confirm whether exhaustive search recovers the simulated occupancy
os.chdir(out_path)
plot_3d_fofc_occ(0.05, 0.95, step=0.05, set_b=40, dataset_prefix=dataset_prefix)


os.chdir(out_path)
for simul_occ in np.arange(0.05, 0.95, 0.05):
    csv_name = "occ_{}_b_{}_u_iso".format(str(simul_occ).replace(".", "_"),set_b)
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
