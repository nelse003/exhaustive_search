# Validation routines that aren't curently used. Will need adaptation to run again.


def occ_loop_merge_confs_simulate_with_refmac_0(bound_state_pdb_path,
                                             ground_state_pdb_path,
                                             input_mtz,
                                             dataset_prefix,
                                             out_path,
                                             set_b = None):

    """ Run exhaustive search on files after refmac 0 cycles"""

    for lig_occupancy in np.arange(0.05, 0.96, 0.05):

        merged_pdb = os.path.join(out_path,"{}_refine_occ_{}.pdb".format(
            dataset_prefix,str(lig_occupancy).replace(".", "_")))

        os.system("giant.merge_conformations input.major={} input.minor={} "
                  "major_occupancy={} minor_occupancy={} output.pdb={}".format(
            ground_state_pdb_path, bound_state_pdb_path, str(1 - lig_occupancy), str(lig_occupancy), merged_pdb))

        simulate_log = os.path.join(out_path,
                                    "{}_simul_{}.log".format(
                                        dataset_prefix,
                                        str(lig_occupancy).replace(".", "_")))

        simulate_mtz = os.path.join(out_path,
                                    "{}_simul_{}.mtz".format(
                                        dataset_prefix,
                                        str(lig_occupancy).replace(".", "_")))

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

            file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive.py"
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
