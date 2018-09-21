from __future__ import division, print_function
import os
import time
import sqlite3
import numpy as np
import pandas as pd
from iotbx.pdb import hierarchy

from select_atoms import get_occupancy_groups

def get_xtals_from_db(params,
                      refinement_outcomes="'3 - In Refinement',"
                                          "'4 - CompChem ready', "
                                          "'5 - Deposition ready',"
                                          "'6 - Deposited'"):

    assert os.path.isfile(params.input.database_path), \
        "The database file: \n {} \n does not exist".format(params.input.database_path)

    # Open connection to sqlite database
    conn = sqlite3.connect(params.input.database_path)
    cur = conn.cursor()

    cur.execute("SELECT CrystalName, RefinementPDB_latest, RefinementMTZ_latest "
                "FROM mainTable WHERE RefinementOutcome in ({})" 
                " AND  (RefinementPDB_latest AND RefinementMTZ_latest) IS NOT NULL".format(refinement_outcomes))

    refinement_xtals = cur.fetchall()

    # Close connection to the database
    cur.close()

    for xtal_name, pdb, mtz in refinement_xtals:
        pdb = pdb.encode('ascii')
        mtz = mtz.encode('ascii')
        xtal_name = xtal_name.encode('ascii')
        yield xtal_name, pdb, mtz

def chunks(l, n):

    """ Divide a list l into chunks of length n. Yield with consecutive letters
    
    Used for splitting a list of atomic points into breaks of 9999 for 
    display in a pdb file. 
    
    Parameters
    -------------------
    l: list
        list to be split into chunks
    y: int
        length to chunk list into
    
    Yields
    ----------------------
    list
        slice of original list up to chunk size
    str
        A letter associated with chunk (A-Z)
    """

    alphabet = []
    for letter in range(65, 91):
        alphabet.append(chr(letter))

    # For item i in a range that is a length of l,
    pos = 0
    for i in xrange(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n], alphabet[pos]
        pos += 1

def write_pdb_HOH_site_cart(pdb, sites_cart):

    """ Write a PDB file containing waters at supplied cartesian sites
     
    Parameters
    -------------------
    pdb: str
        Path to input PDB that contains the crystal record for the 
        crystal symmertry
    sites_cart:
        
    """

    pdb_in = hierarchy.input(file_name=params.input.pdb)

    f = open(params.input.pdb)
    for line in f:
        if line.startswith("CRYST1"):
            cryst = line
        if line.startswith("SCALE1"):
            scale1 = line
        if line.startswith("SCALE2"):
            scale2 = line
        if line.startswith("SCALE3"):
            scale3 = line

    f.close()

    f_out = open("sites_cart.pdb", "w")
    f_out.write(cryst)
    f_out.write(scale1)
    f_out.write(scale2)
    f_out.write(scale3)

    for sites, chain in chunks(sites_cart, 9999):
        print(chain, len(sites))

        for i,site in enumerate(sites):
            f_out.write("HETATM{:>5}  O   HOH {}{:>4}{:>12.3f}{:>8.3f}{:>8.3f}  1.00 10.00           O\n".format(
               i+1, chain, i+1,site[0], site[1], site[2]))

    f_out.close()

def process_validation_csvs(start_occ,
                     end_occ,
                     step,
                     set_b,
                     out_dir,
                     params):

    min_fofcs = []
    min_occs = []
    min_b_facs = []
    fofcs = []
    occs = []
    b_facs = []

    for lig_occupancy in np.arange(start_occ, end_occ + (step / 5), step):
        # TODO Replace CSV naming #59

        csv_name = params.exhaustive.output.csv_prefix \
                   + "_occ_{}_b_{}.csv".format(
            str(lig_occupancy).replace(".", "_"),
            str(set_b).replace(".", "_"))

        csv_path = os.path.join(out_dir, csv_name)
        min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_path)
        fofc = get_fofc_from_csv(csv_path, lig_occupancy,
                                 round_step(b_to_u_iso(set_b)),
                                 step)
        fofcs.append(fofc)
        occs.append(lig_occupancy)
        b_facs.append(set_b)
        min_b_facs.append(u_iso_to_b_fac(min_u_iso))
        min_fofcs.append(fo_fc_at_min)
        min_occs.append(min_occ)

    return min_fofcs, min_occs, min_b_facs, fofcs, occs, b_facs


def expand_array(array):
    "Expand a 3d numpy array"

    x = array[:, 0]
    y = array[:, 1]
    z = array[:, 2]

    return x, y, z


def sample_spherical(npoints, ndim=3):
    """Sample a ndimensional sphere using gaussians
    
    Is used to generate points within a sphere from which 
    a convex hull around atoms can be generated
    """

    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

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
    return round(base * round(float(x)/base), prec)


def set_b_fac_all_occupancy_groups(input_pdb, output_pdb, b_fac, params):

    """ Change b factor of all atoms involved in occupancy groups"""

    pdb_inp = hierarchy.input(input_pdb)
    occ_group = get_occupancy_groups(pdb=input_pdb, params=params)
    for group in occ_group[0]:
        for residue in group:
            for chain in pdb_inp.hierarchy.only_model().chains():
                if chain.id == residue.get('chain'):
                    for residue_group in chain.residue_groups():
                        if residue_group.resseq == residue.get('resseq'):
                            for atom_group in residue_group.atom_groups():
                                for atom in atom_group.atoms():
                                    atom.set_b(b_fac)
                                    print("Changed: {}".format(residue))

    with open(output_pdb, "w") as f:
        f.write(pdb_inp.hierarchy.as_pdb_string(
            crystal_symmetry=pdb_inp.input.crystal_symmetry()))


def set_b_fac_all_atoms(input_pdb, output_pdb, b_fac):

    """ Change B factors of all atoms to the same provided value"""
    pdb_inp = hierarchy.input(input_pdb)
    for chain in pdb_inp.hierarchy.only_model().chains():
        for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                for atom in atom_group.atoms():
                    atom.set_b(b_fac)

    with open(output_pdb, "w") as f:
        f.write(pdb_inp.hierarchy.as_pdb_string(
            crystal_symmetry=pdb_inp.input.crystal_symmetry()))


def get_minimum_fofc(csv_name, b_fac=None):
    """
    Get minima in fofc, and return minima and where it occurs

    B factor can be supplied to look at the minima across a
    single b factor value

    :param csv_name:
    :param b_fac:
    :return:
    """
    print(os.getcwd())
    # TODO Remove this dual .csv by cleaning up csv name: issue 59
    if csv_name.endswith(".csv"):
        data = np.genfromtxt(csv_name, delimiter=',', skip_header=0)
    else:
        data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',',
                             skip_header=0)

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
        raise ValueError("|Fo-Fc| data is not in correct format")

    if b_fac is not None:
        set_u_iso = b_to_u_iso(b_fac)
        step = np.unique(occ)[1] - np.unique(occ)[0]
        check_u_iso = round_step(set_u_iso, base=step)
        data_array = np.stack((occ, u_iso, fo_fc))

        occ = data_array[0][(data_array[1] >= check_u_iso)
                            & (data_array[1] <= check_u_iso)]

        u_iso = data_array[1][(data_array[1] >= check_u_iso)
                              & (data_array[1] <= check_u_iso)]

        fo_fc = data_array[2][(data_array[1] >= check_u_iso)
                              & (data_array[1] <= check_u_iso)]

    min_index = np.argmin(fo_fc)

    return occ[min_index], u_iso[min_index], fo_fc[min_index]


def get_lig_occ(refine_pdb):

    pdb_in = hierarchy.input(refine_pdb)

    lig_atoms = []
    for chain in pdb_in.hierarchy.only_model().chains():
        for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                if atom_group.resname == "LIG":
                    for atom in atom_group.atoms():
                        lig_atoms.append((atom_group.altloc, atom.occ))

    if len(list(set(lig_atoms))) == 2:
        end_occ = list(set(lig_atoms))[0][1]+list(set(lig_atoms))[1][1]
        return end_occ
    else:
        "Ligand occupancy is not defined"
        exit()


def wait_for_file_existence(file_path, wait_time):

    """ Wait for a file to exist, stop after wait_time (seconds)

    Used for waiting for qsub to finish
    """

    time_in_loop = 0
    while not os.path.exists(file_path):
        if time_in_loop < wait_time:
            print("waiting")
            time.sleep(1)
            time_in_loop += 1
        else:
            raise IOError("Cannot find file {} within {} seconds".format(
                file_path, wait_time))


def get_csv_filepath(params):
    for occupancy in np.arange(params.validate.options.start_simul_occ,
                               params.validate.options.end_simul_occ +
                               params.validate.options.step_simulation / 5,
                               params.validate.options.step_simulation):

        yield os.path.join(params.output.out_dir,
                           params.exhaustive.output.csv_prefix +
                           "_occ_{}_b_{}.csv".format(
                               str(occupancy).replace(".", "_"),
                               str(params.validate.options.set_b).replace(".", "_")))


def get_random_starting_occ_from_folder_name(occupancy, out_path,
                                             dataset_prefix):

    """Pull occupancy values from folder names of random refinements"""

    folders = [name for name in os.listdir(
        os.path.join(out_path, dataset_prefix
                     + "_refine_occ_"
                     + str(occupancy).replace(".", "_"))) if
               os.path.isdir(
                   os.path.join(out_path, dataset_prefix
                                + "_refine_occ_"
                                + str(occupancy).replace(".", "_"),
                                name))]

    folders = [name for name in folders if "exhaustive" not in name]

    for folder in folders:
        yield float("0." + folder.split('_')[-1])


def read_ligand_occupancy_b(pdb_path, lig_chain):
    """Extract occupancy and B factor of ligand of
    interest from one PDB file into a dataframe"""

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

    occ_b_df = pd.DataFrame(lig_occ_b,
                            columns=["Atom", "Occupancy", "B_factor"])

    return occ_b_df


def get_pdbs(refinement_dir, pdb_name="refine.pdb"):
    """ Given a folder get all pdb paths that match name"""

    pdbs = []
    for root, dirs, files in os.walk(refinement_dir):
        for f in files:
            if f == pdb_name:
                pdbs.append(os.path.join(root, f))

    print("# of pdbs", len(pdbs))
    return pdbs


def get_occ_b(refinement_dir, lig_chain,
              pdb_name="refine.split.bound_state.pdb"):
    """ Given a folder get occupancies & B factor of ligand in chain"""

    pdbs = get_pdbs(refinement_dir, pdb_name)

    occ_b = []

    for pdb in pdbs:

        dataset, _ = os.path.split(pdb)
        dataset = dataset.split("/")[-1]

        occ_b_df = read_ligand_occupancy_b(pdb, lig_chain)
        mean_ligand_b_factor = occ_b_df['B_factor'].mean()
        std_ligand_b_factor = occ_b_df['B_factor'].std()

        if occ_b_df.apply(lambda x: x.nunique())[1] == 1:
            lig_occ = occ_b_df.loc("Occupancy")[0][1]
            occ_b.append([dataset, lig_occ,
                          mean_ligand_b_factor,
                          std_ligand_b_factor])
        else:
            print(occ_b_df)
            print("Occupancy varies across ligand, "
                  "histogram not currently generated")

    print(occ_b)
    occ_df = pd.DataFrame(data=occ_b,
                          columns=['dataset', 'occupancy', 'mean_b_fac',
                                   'std_b_fac'])

    return occ_df


def get_fofc_from_csv(csv_name, occupancy, u_iso, step=0.05):

    occupancy = round_step(occupancy, base=step)
    u_iso = round_step(u_iso, base=step)

    #TODO Remove this dual .csv by cleaning up csv name: issue 59

    if csv_name.endswith(".csv"):
        data = np.genfromtxt(csv_name, delimiter=',', skip_header=0)
    else:
        data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',',
                             skip_header=0)

    e=0.0001
    data_line = data[((occupancy-e)<data[:,0])
                     & ((occupancy+e)>data[:,0])
                     & ((u_iso-e)<data[:,2])
                     & ((u_iso+e)>data[:,2]) ]

    print(u_iso)
    print(data_line)

    fo_fc = data_line[0][3]

    return fo_fc


def datasets_from_compound(protein_prefix,compound_folder):

    """ Loop over datasets with prefix in folder."""

    for dataset in filter(lambda x:
                          x.startswith(protein_prefix)
                          and os.path.isdir(os.path.join(compound_folder, x)),
                                os.listdir(compound_folder)):
        yield dataset


def collate_edstats_scores(protein_prefix, compound_folder):

    """ Collate edstats scores into DataFrame. Write csv. Return df"""

    dfs = []
    for dataset in datasets_from_compound(protein_prefix, compound_folder):

        edstats_csv = os.path.join(compound_folder, dataset, "Plots",
                                   "residue_scores.csv")

        if os.path.exists(edstats_csv):
            edstats_df = pd.read_csv(edstats_csv)
            edstats_df['Dataset'] = dataset
            dfs.append(edstats_df)
        else:
            print("No edstats scores found for {}".format(dataset))
    try:
        compound_edstats = pd.concat(dfs, ignore_index=True)
    except ValueError:
        print("Edstats not able to concatanate")
        return None

    print(compound_edstats)

    compound_edstats.to_csv(file_name=os.path.join(compound_folder,
                                                   "collated_edstats"),
                            index=False)

    return compound_edstats


def remove_residues(input_pdb, output_pdb, residues_remove):

    pdb_in = hierarchy.input(file_name = input_pdb)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()

    selection_string_list = []
    for residue_remove in residues_remove:
        selection_string = "(resid {} and chain {} and altid {})".format(residue_remove[0],
                                                            residue_remove[1],residue_remove[2])
        selection_string_list.append(selection_string)

    selection_string = " or ".join(selection_string_list)
    not_selection_string ="not ({})".format(selection_string)

    acceptor_hierarchy = pdb_in.construct_hierarchy()

    print(not_selection_string)

    remove_atoms_sel = sel_cache.selection(not_selection_string)
    removed_hier = acceptor_hierarchy.select(remove_atoms_sel)

    f = open(os.path.join(output_pdb), "w+")

    f.write(removed_hier.as_pdb_string(
        crystal_symmetry=pdb_in.input.crystal_symmetry()))

    f.close()


def is_almost_equal(x,y, epsilon=1*10**(-8)):
    """Return True if two values are close in numeric value
        By default close is withing 1*10^-8 of each other
        i.e. 0.00000001
    """
    return abs(x-y) <= epsilon
