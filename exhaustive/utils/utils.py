import os
import time
import numpy as np
import pandas as pd
from iotbx.pdb import hierarchy

from exhaustive.utils.select import get_occupancy_groups

def b_to_u_iso(b_fac):
    """ Convert isotropic B factor to u iso"""

    u_iso = b_fac/(8 * np.pi ** 2)
    return u_iso

def u_iso_to_b_fac(u_iso):
    """ Convert u_iso to isotropic B factor """

    b_iso = (8 * np.pi ** 2) * u_iso
    return b_iso

def round_step(x, prec=2, base=.05):
    """ Return a number rounded to the nearest base."""
    return round(base * round(float(x)/base),prec)

def get_fofc_from_csv(csv_name,occupancy, u_iso, step=0.05):

    occupancy = round_step(occupancy,base=step)
    u_iso = round_step(u_iso,base=step)

    #TODO Remove this dual .csv by cleaning up csv name #59

    if csv_name.endswith(".csv"):
        data = np.genfromtxt(csv_name, delimiter=',', skip_header=0)
    else:
        data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)

    e=0.0001
    data_line = data[((occupancy-e)<data[:,0]) & ((occupancy+e)>data[:,0]) & ((u_iso-e)<data[:,2]) & ((u_iso+e)>data[:,2]) ]

    print u_iso
    print data_line

    fo_fc = data_line[0][3]

    return fo_fc

def set_b_fac_all_occupancy_groups(input_pdb, output_pdb, b_fac, params):

    """ Change b factor of all atoms involved in occupancy groups"""

    pdb_inp = hierarchy.input(input_pdb)
    occ_group = get_occupancy_groups(pdb=input_pdb, params = params)
    for group in occ_group[0]:
        for residue in group:
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

def set_b_fac_all_atoms(input_pdb, output_pdb, b_fac):

    """ Change B factors of all atoms to the same provided value"""
    pdb_inp = hierarchy.input(input_pdb)
    for chain in pdb_inp.hierarchy.only_model().chains():
        for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                for atom in atom_group.atoms() :
                    atom.set_b(b_fac)
    with open(output_pdb, "w") as f:
        f.write(pdb_inp.hierarchy.as_pdb_string(crystal_symmetry=pdb_inp.input.crystal_symmetry()))

def get_minimum_fofc(csv_name, b_fac=None):
    """
    Get minima in fofc, and return minima and where it occurs

    B factor can be supplied to look at the minima across a single b factor value

    :param csv_name: 
    :param b_fac: 
    :return: 
    """
    print(os.getcwd())
    #TODO Remove this dual .csv by cleaning up csv name #59
    if csv_name.endswith(".csv"):
        data = np.genfromtxt(csv_name, delimiter=',', skip_header=0)
    else:
        data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)

    # If four column data from multiple ligand
    if len(data[0]) == 4:
        occ = data[:, 0]
        u_iso = data[:, 2]
        fo_fc = data[:, 3]
    elif len(data[0]) == 3:
        occ = data[:, 0]
        u_iso = data[:, 1]
        fo_fc = data[:, 2]
    else:
        print("Data is not in correct format")
        exit()

    if b_fac is not None:
        set_u_iso = b_to_u_iso(b_fac)
        step = np.unique(occ)[1] - np.unique(occ)[0]
        check_u_iso = round_step(set_u_iso, base=step)
        data_array = np.stack((occ, u_iso, fo_fc))
        occ = data_array[0][(data_array[1] >= check_u_iso) & (data_array[1] <= check_u_iso)]
        u_iso = data_array[1][(data_array[1] >= check_u_iso) & (data_array[1] <= check_u_iso)]
        fo_fc = data_array[2][(data_array[1] >= check_u_iso) & (data_array[1] <= check_u_iso)]

    min_index = np.argmin(fo_fc)

    return occ[min_index], u_iso[min_index], fo_fc[min_index]

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

def wait_for_file_existence(file_path, wait_time):

    """ Wait for a file to exist, stop after wait_time (seconds)

    Used for waitign for qsub to finish
    """

    time_in_loop = 0
    while not os.path.exists(file_path):
        if time_in_loop < wait_time:
            print("waiting")
            time.sleep(1)
            time_in_loop += 1
        else:
            raise IOError("Cannot find file {} within {} seconds".format(file_path, wait_time))

def get_csv_filepath(params):
    for occupancy in np.arange(params.validate.options.start_simul_occ,
                               params.validate.options.end_simul_occ +
                               params.validate.options.step_simulation / 5,
                               params.validate.options.step_simulation):

        yield os.path.join(params.output.out_dir, params.exhaustive.output.csv_prefix + \
                                            "_occ_{}_b_{}.csv".format(str(occupancy).replace(".", "_"),
                                                                      str(params.validate.options.set_b).replace(".","_")))

def get_random_starting_occ_from_folder_name(occupancy, out_path, dataset_prefix):

    "From folder structure pull out occupancy values from random refinements"

    folders = [name for name in os.listdir(
        os.path.join(out_path, dataset_prefix + "_refine_occ_" + str(occupancy).replace(".", "_"))) if
               os.path.isdir(
                   os.path.join(out_path, dataset_prefix + "_refine_occ_" + str(occupancy).replace(".", "_"),
                                name))]

    folders = [name for name in folders if "exhaustive" not in name]

    for folder in folders:
        yield float("0." + folder.split('_')[-1])

###############################################################################
#  Likely depreceated.
###############################################################################

def get_delta_fofc(csv_name, occupancy, u_iso, step=0.05):

    """ Find the difference between mean(|Fo-Fc|) at minima and at supplied u_iso and occ"""

    fo_fc = get_fofc_from_csv(csv_name,occupancy, u_iso)

    min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_name)
    delta_fofc = fo_fc - fo_fc_at_min

    return delta_fofc

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

def read_ligand_occupancy_b(pdb_path, lig_chain):
    """Extract occupancy and B factor of ligand of interest from one PDB file into a dataframe"""

    # Read in single PDB file
    pdb_in = hierarchy.input(file_name=pdb_path)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()
    lig_sel = sel_cache.selection("chain {}".format(lig_chain))
    lig_hierarchy = pdb_in.hierarchy.select(lig_sel)

    lig_occ_b = []
    # Get occupancy & B factor of ligand
    for model in lig_hierarchy.models():
        for chain in model.chains():
            for rg in chain.residue_groups():
                for ag in rg.atom_groups():
                    for atom in ag.atoms():
                        lig_occ_b.append([atom.name, atom.occ, atom.b])
    occ_b_df = pd.DataFrame(lig_occ_b, columns=["Atom", "Occupancy", "B_factor"])

    return occ_b_df

def get_pdbs(refinement_dir, pdb_name="refine.pdb"):
    """ Given a folder get all pdb paths that match compound and/or name

    If a compound code is provided will only extract pdbs with that compound cif present in folder
    """

    pdbs = []
    for root, dirs, files in os.walk(refinement_dir):
        for file in files:
            if file == pdb_name:
                pdbs.append(os.path.join(root, file))


    print("# of pdbs", len(pdbs))
    return pdbs


def get_occs(refinement_dir, lig_chain, pdb_name="refine.split.bound_state.pdb"):
    """ Given a folder get all occupancies of ligand in provided chain

    If a compound code is provided will only extract pdbs with that compound cif present in folder
    """

    pdbs = get_pdbs(refinement_dir, pdb_name)

    all_lig_occupancy = []
    for pdb in pdbs:

        occ_b_df = read_ligand_occupancy_b(pdb, lig_chain)

        if occ_b_df.apply(lambda x: x.nunique())[1] == 1:
            lig_occ = occ_b_df.loc("Occupancy")[0][1]
            all_lig_occupancy.append(lig_occ)
        else:
            print(occ_b_df)
            print "Occupancy varies across ligand, histogram not currently generated"

    return all_lig_occupancy