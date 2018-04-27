import os
import random
import time

import iotbx.mtz
import matplotlib.pyplot as plt
import numpy as np
from iotbx.pdb import hierarchy

from plot_exhaustive_search import scatter_plot
from process_exhaustive_search import get_minimum_fofc
from refinement import refmac_0_cyc
from select_occupancy_groups import get_occupancy_groups
from utils import b_to_u_iso, round_step


def occ_loop_merge_confs_simulate(bound_state_pdb_path,
                                 ground_state_pdb_path,
                                 input_mtz,
                                 dataset_prefix,
                                 out_path,
                                 set_b = None,
                                 step = 0.05,
                                 start_occ = 0.05,
                                 end_occ = 0.95,
                                 buffer = 0,
                                 grid_spacing = 0.25):

    for lig_occupancy in np.arange(start_occ, end_occ+step/5, step):
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
        csv_name = "occ_{}_b_{}_u_iso".format(str(lig_occupancy).replace(".", "_"),str(set_b).replace(".", "_"))

        with open(os.path.join(out_path, sh_file),'w') as file:

            file.write("#!/bin/bash\n")
            file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
            file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")

            file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive_search.py"
                       " input.pdb={} input.mtz={} output.out_dir={} xtal_name={} "
                       "options.csv_name={} options.step={} options.buffer={} "
                       "options.grid_spacing={}".format(merged_pdb, merged_pdb+".mtz",out_path, dataset_prefix,csv_name,
                                                         step, buffer, grid_spacing))

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

def submit_exhasutive_with_refmac_0(dataset_prefix, out_path, set_b = None):

    """ Qsub submission of sh files from occ_loop_merge_confs_simulate_with_refmac_0"""

    for lig_occupancy in np.arange(0.05, 0.96, 0.05):

        sh_file = "{}_occ_{}_b_{}.sh".format(dataset_prefix,
                                             str(lig_occupancy).replace(".", "_"),
                                             str(set_b).replace(".", "_"))

        os.system("qsub -o {} -e {} {}".format(
            os.path.join(out_path, "output_{}.txt".format(str(lig_occupancy).replace(".", "_"))),
            os.path.join(out_path, "error_{}.txt".format(str(lig_occupancy).replace(".", "_"))),
            os.path.join(out_path, sh_file)))

def set_u_iso_all_occupancy_groups(input_pdb, output_pdb, b_fac):

    pdb_inp = hierarchy.input(input_pdb)
    occ_group = get_occupancy_groups(pdb= input_pdb)
    for group in occ_group[0]:
        for residue in group:

            #TODO Replace many for loops with better structure

            for chain in pdb_inp.hierarchy.only_model().chains():
                if chain.id == residue.get('chain'):
                    for residue_group in chain.residue_groups():
                        if residue_group.resseq == residue.get('resseq'):
                            for atom_group in residue_group.atom_groups():
                                for atom in atom_group.atoms() :
                                    atom.set_b(b_fac)
                                    print("Changed: {}".format(residue))
    with open(output_pdb,"w") as f:
        f.write(pdb_inp.hierarchy.as_pdb_string(crystal_symmetry=pdb_inp.input.crystal_symmetry()))

def read_occ_csv(csv_name):

    """Return numpy arrays from  a csv with lig_occ, ground_occ, u_iso, mean(|Fo-Fc|)"""

    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)
    occ = data[:, 0]
    u_iso = data[:, 1]
    fo_fc = data[:, 2]

    return occ, u_iso, fo_fc

def get_delta_occ(csv_name, occupancy):

    """ Find the difference between minima occupancy and supplied occupancy.
     
     Given a csv with lig_occ, ground_occ, u_iso, mean(|Fo-Fc|) csv and occupancy value.
     """

    min_occ, _, _ = get_minimum_fofc(csv_name)
    return occupancy - min_occ

def get_delta_u_iso(csv_name):

    """ Find the minima u_iso and supplied u_iso.

     Given a csv with lig_occ, ground_occ, u_iso, mean(|Fo-Fc|) csv and u_iso value. 
     Useful when using a set B factor.
     """

    _, min_u_iso, _ = get_minimum_fofc(csv_name)
    return u_iso - min_u_iso

def get_delta_fofc(csv_name, occupancy, u_iso, step=0.05):

    """ Find the difference between mean(|Fo-Fc|) at minima and at supplied u_iso and occ"""

    fo_fc = get_fofc_from_csv(csv_name,occupancy, u_iso)

    min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_name)
    delta_fofc = fo_fc - fo_fc_at_min

    return delta_fofc

def get_fofc_from_csv(csv_name,occupancy, u_iso, step=0.05):

    occupancy = round_step(occupancy,base=step)
    u_iso = round_step(u_iso,base=step)

    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)
    e=0.0001
    data_line = data[((occupancy-e)<data[:,0]) & ((occupancy+e)>data[:,0]) & ((u_iso-e)<data[:,2]) & ((u_iso+e)>data[:,2]) ]
    fo_fc = data_line[0][3]

    return fo_fc

def get_delta_fofc_over_occupancies(set_b,start_occ, end_occ, step = 0.05):

    """ Get delta occupancy across occupancies"""

    occupancies = []
    delta_fofcs = []

    for lig_occupancy in np.arange(start_occ, end_occ+(step/5), step):
        delta_fofc = get_delta_fofc("occ_{}_b_{}_u_iso".format(
                                    str(lig_occupancy).replace(".","_"),
                                    str(set_b).replace(".","_")),
                                    lig_occupancy,
                                    b_to_u_iso(set_b),
                                    step = step)
        occupancies.append(lig_occupancy)
        delta_fofcs.append(delta_fofc)

    return np.column_stack((occupancies,delta_fofcs))

def connectpoints(x,y,x_1,y_1,p1):
    x1, x2 = x[p1], x_1[p1]
    y1, y2 = y[p1], y_1[p1]
    plt.plot([x1,x2],[y1,y2],'k-')

def connectpoint(x,y,x_1,y_1):
    plt.plot([x,x_1],[y,y_1],'k-')

def plot_fofc_occ(start_occ, end_occ, step, set_b):

    min_fofcs = []
    min_occs = []
    fofcs = []
    occs = []

    for lig_occupancy in np.arange(start_occ, end_occ + (step / 5), step):
        csv_name = "occ_{}_b_{}_u_iso".format(str(lig_occupancy).replace(".","_"),str(set_b).replace(".","_"))
        min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_name)
        fofc = get_fofc_from_csv(csv_name,lig_occupancy, round_step(b_to_u_iso(40)), step)
        fofcs.append(fofc)
        occs.append(lig_occupancy)
        min_fofcs.append(fo_fc_at_min)
        min_occs.append(min_occ)

    fig, ax = plt.subplots()
    min_plot, = ax.plot(min_occs, min_fofcs,'k+')
    occ_plot, = ax.plot(occs, fofcs, 'ro')

    for i in np.arange(0, len(occs)):
        connectpoints(occs,fofcs,min_occs,min_fofcs,i)

    ax.legend((min_plot,occ_plot),
              ('Minima of mean |Fo-Fc|','Mean |Fo-Fc| at simulated occupancy'),
              prop={"size": 8},
              numpoints=1,
              bbox_to_anchor=(1, 1),
              bbox_transform=plt.gcf().transFigure)

    ax.set_xlabel("Occupancy")
    ax.set_ylabel("Mean |Fo-Fc|")

    plt.title("Delta mean|Fo-Fc| and Delta Occupancy", fontsize=10)
    plt.savefig("NUDT7A-x1740-set-b-delta_fofc_occ.png")

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

def quick_refine_repeats(start_occ, end_occ,step, dataset_prefix, set_b, out_path, input_cif):

    params = "multi-state-restraints.refmac.params"

    for simul_occ in np.arange(start_occ, end_occ + step / 5, step):
        for starting_rand_occ in get_starting_occ(simul_occ, out_path, dataset_prefix):

            simulate_mtz = os.path.join(out_path,"{}_refine_occ_{}".format(dataset_prefix,str(simul_occ).replace(".", "_")),
                                        "{}_simul_{}.mtz".format(dataset_prefix, str(simul_occ).replace(".", "_")))

            out_path = os.path.join(out_path,"{}_refine_occ_{}".format(dataset_prefix,str(simul_occ).replace(".", "_")),
                                    "{}_expected_occ_{}_b_{}_supplied_occ_{}".format(dataset_prefix,
                                                                                     str(simul_occ).replace(".", "_"),
                                                                                     str(set_b).replace(".","_"),
                                                                                     str(starting_rand_occ).replace(".", "_")))
            if not os.path.exists(out_path):
                os.mkdir(out_path)
            os.chdir(out_path)

            refinement_random_occ_pdb = os.path.join(out_path,"{}_random_refine_occ_{}.pdb".format(
                                                                    dataset_prefix,
                                                                    str(starting_rand_occ).replace(".", "_")))

            sh_file = "{}_expected_occ_{}_b_{}_supplied_occ_{}.sh".format(dataset_prefix,
                                                                          str(simul_occ).replace(".", "_"),
                                                                          str(set_b).replace(".", "_"),
                                                                          str(starting_rand_occ).replace(".", "_"))

            out_prefix = "expected_occ_{}_supplied_occ_{}".format(str(simul_occ).replace(".", "_"),
                                                                  str(starting_rand_occ).replace(".", "_"))
            dir_prefix = "refine_" + out_prefix + "_cyc_20"

            print(out_path)

            os.system("cp multi-state-restraints.refmac.params multi-state-restraints-tmp.refmac.params")
            with open("multi-state-restraints-tmp.refmac.params",'a') as file:
                file.write('ncyc 20')

            with open(os.path.join(out_path, sh_file), 'w') as file:
                file.write("#!/bin/bash\n")
                file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
                file.write(
                    "source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
                file.write("giant.quick_refine input.pdb={} input.mtz={} input.cif={} "
                           "input.params={} output.out_prefix={} output.dir_prefix={} ".format(
                    refinement_random_occ_pdb, simulate_mtz, input_cif, "multi-state-restraints-tmp.refmac.params",
                    out_prefix, dir_prefix,))

            os.system("qsub -o {} -e {} {}".format(
                os.path.join(out_path, "output_{}.txt".format(str(simul_occ).replace(".", "_"))),
                os.path.join(out_path, "error_{}.txt".format(str(simul_occ).replace(".", "_"))),
                os.path.join(out_path, sh_file)))

            out_path = os.path.dirname(os.path.dirname(out_path))

def refine_after_exhasutive_search(input_pdb, input_mtz, input_cif, refine_params, dataset_prefix, working_dir = None):

    os.chdir(working_dir)
    sh_file = "{}_quick_refine_exhaustive_search_minima.sh".format(dataset_prefix)
    out_prefix = "{}_refine_after_exhaustive_search".format(dataset_prefix)
    dir_prefix = out_prefix
    quick_refine_qsub(input_pdb=input_pdb, input_mtz=input_mtz, input_cif=input_cif, refine_params=refine_params,
                      sh_file=sh_file, out_prefix=out_prefix, dir_prefix=dir_prefix)


def quick_refine_qsub(input_pdb, input_mtz, input_cif, refine_params,
                      sh_file, out_prefix=None, dir_prefix=None):

    with open(os.path.join(out_path, sh_file), 'w') as file:
        file.write("#!/bin/bash\n")
        file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
        file.write(
            "source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
        file.write("giant.quick_refine input.pdb={} input.mtz={} input.cif={} "
                   "input.params={} output.out_prefix={} output.dir_prefix={} ".format(input_pdb, input_mtz,
                                                                                       input_cif, refine_params,
                                                                                       out_prefix, dir_prefix))

    os.system("qsub -o {} -e {} {}".format(
            os.path.join(out_path, "output_{}.txt".format(str(out_prefix))),
            os.path.join(out_path, "error_{}.txt".format(str(out_prefix))),
            os.path.join(out_path, sh_file)))

def get_starting_occ(occupancy,out_path, dataset_prefix):

    folders = [name for name in os.listdir(
        os.path.join(out_path, dataset_prefix + "_refine_occ_" + str(occupancy).replace(".","_"))) if
        os.path.isdir(os.path.join(out_path, dataset_prefix + "_refine_occ_" + str(occupancy).replace(".","_"),name))]

    folders = [name for name in folders if "exhaustive" not in name]

    for folder in folders:
       yield float("0." + folder.split('_')[-1])

def get_lig_occ(refine_pdb):

    pdb_in = hierarchy.input(refine_pdb)

    lig_atoms =[]

    for chain in pdb_in.hierarchy.only_model().chains() :
        for residue_group in chain.residue_groups() :
            for atom_group in residue_group.atom_groups() :
                if atom_group.resname == "LIG":
                    for atom in atom_group.atoms() :
                        lig_atoms.append((atom_group.altloc, atom.occ))
    if len(list(set(lig_atoms))) == 2:
        end_occ = list(set(lig_atoms))[0][1]+list(set(lig_atoms))[1][1]
        return end_occ
    else:
        "Ligand occupancy is not defined"
        exit()


def plot_random_refinement_with_ES(start_occ, end_occ,step, dataset_prefix, set_b,out_path):

    for simul_occ in np.arange(start_occ,end_occ+step/5,step):

        refine_5_cyc_occs = []
        refine_20_cyc_occs = []
        rand_occs = []

        # Get occupancy from exhaustive search and refined exhaustive search pdbs

        ES_folder = os.path.join(out_path,dataset_prefix + "_refine_occ_" + str(simul_occ).replace(".", "_"))
        ES_pdb = os.path.join(ES_folder,"exhaustive_search_minima.pdb")

        folders = [name for name in os.listdir(ES_folder) if os.path.isdir(os.path.join(ES_folder, name))]

        es_refine = []
        for folder in folders:
            if folder.find('refine_after_exhaustive_search') != -1:
                es_refine.append(int(folder[-4:]))
        es_refine_folder = os.path.join(ES_folder,
                                       "{}_refine_after_exhaustive_search{}".format(dataset_prefix,
                                                                                    str(max(es_refine)).rjust(4,'0')))
        ES_refine_pdb = os.path.join(es_refine_folder, "{}_refine_after_exhaustive_search.pdb".format(dataset_prefix))

        ES_occ = get_lig_occ(ES_pdb)
        ES_refine_occ = get_lig_occ(ES_refine_pdb)

        # Loop over starting occupancy, got from folder name
        # Get refined occupancy from pdb file
        for starting_rand_occ in get_starting_occ(simul_occ,out_path,dataset_prefix):
            cur_dir = os.path.join(out_path,
                                    dataset_prefix + "_refine_occ_" + str(simul_occ).replace(".", "_"),
                                    dataset_prefix + "_expected_occ_"
                                    + str(simul_occ).replace(".", "_") + "_b_" + str(set_b) + "_supplied_occ_" +
                                    str(starting_rand_occ).replace(".", "_"))

            rand_occs.append(starting_rand_occ)
            #refine_pdb = os.path.join(cur_dir,"refine.pdb")



            # get folders for each refinement cycle
            out_prefix = "refine_expected_occ_{}_supplied_occ_{}".format(str(simul_occ).replace(".", "_"),
                                                                  str(starting_rand_occ).replace(".", "_"))

            folders = [name for name in os.listdir(cur_dir) if os.path.isdir(os.path.join(cur_dir, name))]

            cyc_20 = []
            cyc_5 = []
            for folder in folders:
                if folder.find('cyc_20') != -1:
                    cyc_20.append(int(folder[-4:]))
                elif folder.find('cyc') == -1:
                    cyc_5.append(int(folder[-4:]))

            cur_folder_5_cyc = os.path.join(cur_dir, out_prefix + "_" + str(max(cyc_5)).rjust(4,'0'))
            cur_folder_20_cyc = os.path.join(cur_dir, out_prefix + "_cyc_20" + str(max(cyc_20)).rjust(4,'0'))

            for file in os.listdir(cur_folder_5_cyc):
                if file.find('expected') != -1 & file.find('pdb')!=-1:
                    cyc_5_refine_pdb = os.path.join(cur_folder_5_cyc, file)
            for file in os.listdir(cur_folder_20_cyc):
                if file.find('expected') != -1 & file.find('pdb')!=-1:
                    cyc_20_refine_pdb = os.path.join(cur_folder_20_cyc, file)

            refine_occ_cyc5 = get_lig_occ(cyc_5_refine_pdb)
            refine_5_cyc_occs.append(refine_occ_cyc5)

            refine_occ_cyc20 = get_lig_occ(cyc_20_refine_pdb)
            refine_20_cyc_occs.append(refine_occ_cyc20)

        # Plot the random starting occupancy, refined occupancy and simulated occupancy (target)
        fig, ax = plt.subplots()
        #fig.subplots_adjust(right=2.0)


        rand_occ_plot, = ax.plot(rand_occs,'ko')
        refine_5_cyc_plot, = ax.plot(refine_5_cyc_occs, 'ro')
        refine_20_cyc_plot, = ax.plot(refine_20_cyc_occs, 'mo')
        es_plot, = ax.plot(len(rand_occs),ES_occ,'go')
        es_refine_plot, = ax.plot(len(rand_occs),ES_refine_occ,'yo')
        simul_occ_plot, = ax.plot(np.arange(0, len(rand_occs)+1),[simul_occ]*(len(rand_occs)+1), 'bo')

        # Lines to connect random starting occupancy and refined occupancy
        for i in np.arange(0, len(rand_occs)):
             connectpoints(np.arange(0, len(rand_occs)),rand_occs, np.arange(0, len(rand_occs)),refine_5_cyc_occs, i)

        # Lines to connect refined occupancy to simualted occupancy
        for i in np.arange(0, len(rand_occs)):
             connectpoints(np.arange(0, len(rand_occs)),refine_5_cyc_occs,
                           np.arange(0, len(rand_occs)),[simul_occ]*len(rand_occs), i)

        # Lines to connect refined occupancy to simualted occupancy
        for i in np.arange(0, len(rand_occs)):
             connectpoints(np.arange(0, len(rand_occs)),refine_20_cyc_occs,
                           np.arange(0, len(rand_occs)),[simul_occ]*len(rand_occs), i)

        # Line to conenct ES occ and ES refined occ

        connectpoint(len(rand_occs),ES_occ,
                      len(rand_occs),ES_refine_occ)

        # Line to connect ES occ and simulated occ

        connectpoint(len(rand_occs),ES_occ,
                      len(rand_occs),simul_occ)

        lgd = ax.legend((rand_occ_plot, refine_5_cyc_plot, refine_20_cyc_plot, simul_occ_plot,es_plot,es_refine_plot),
                  ('Starting occupancy (random)', 'Refined occupancy (5 cycles)',
                   'Refined Occupancy (20 cycles)','Simulated Occupancy (target)',
                   'Minima of exhaustive search', 'Minima of exhaustive search, refined (5 cycles)' ),
                  prop={"size": 8},
                  numpoints=1,
                  bbox_to_anchor=(1.35, 0.5),
                  bbox_transform=plt.gcf().transFigure)

        ax.set_ylabel("Occupancy")

        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom='off',  # ticks along the bottom edge are off
            top='off',  # ticks along the top edge are off
            labelbottom='off')  # labels along the bottom edge are off

        ax.set_xlabel("Random refinement start points")

        plt.xlim(-1, len(rand_occs) +2)

        plt.savefig(os.path.join(ES_folder,"random_occ_start_refine_simul_{}_20_cyc.png".format(
            str(simul_occ).replace('.','_'))), bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)


def wait_for_file_existence(file_path, wait_time):
    time_in_loop = 0
    while not os.path.exists(file_path):
        if time_in_loop < wait_time:
            print("waiting")
            time.sleep(1)
            time_in_loop += 1
        else:
            raise IOError("Cannot find file {} within {} seconds".format(file_path, wait_time))


def get_csv_filepath(directory, set_b=None, step=0.05, start_occ=0.05, end_occ=0.95):
    for occupancy in np.arange(start_occ, end_occ + step / 5, step):
        if set_b is not None:
            yield os.path.join(directory, "occ_{}_b_{}_u_iso.csv".format(str(occupancy).replace(".","_"), set_b))
        else:
            yield os.path.join(directory, "occ_{}_u_iso.csv".format((str(occupancy).replace(".","_"))))


in_path = "/dls/labxchem/data/2018/lb18145-55/processing/analysis/initial_model/NUDT22A-x1038"
bound_state_pdb_path = os.path.join(in_path, "refine.output.bound-state.pdb")
ground_state_pdb_path = os.path.join(in_path, "refine.output.ground-state.pdb")
input_mtz = os.path.join(in_path, "NUDT22A-x1038.free.mtz")
# input_cif = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1740/NUOOA000181a.cif"
dataset_prefix = "NUDT22A-x1038"
out_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/exhaustive_search_phenix_fmodel/NUDT22A-x1038"

if not os.path.exists(out_path):
    os.mkdir(out_path)

# This loop runs exhaustive search many times across simulated data
occ_loop_merge_confs_simulate(bound_state_pdb_path,
                              ground_state_pdb_path,
                              input_mtz,
                              dataset_prefix,
                              out_path,
                              set_b = 40,
                              step = 0.01,
                              start_occ = 0.01,
                              end_occ = 0.99,
                              buffer = 0,
                              grid_spacing = 0.25)


# Waits for occupancy csvs to be output
for file_path in get_csv_filepath(out_path, set_b=40, step=0.01, start_occ=0.01, end_occ=0.99):
    wait_for_file_existence(file_path, wait_time=1000)

# This plots exhaustive search results, to confirm whether exhaustive search recovers the simulated occupancy
os.chdir(out_path)
plot_fofc_occ(0.01, 0.99, step=0.01, set_b=40)


os.chdir(out_path)
for simul_occ in np.arange(0.01, 0.99, 0.01):
    csv_name = "occ_{}_b_40_u_iso".format(str(simul_occ).replace(".", "_"))
    scatter_plot(csv_name, title_text="Phenix.fmodel at occ {}".format(simul_occ))
exit()
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

    for starting_rand_occ in get_starting_occ(simul_occ, validation_path, dataset_prefix):
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
