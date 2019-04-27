from __future__ import division
from __future__ import print_function

import os
import sqlite3

import numpy as np
import pandas as pd
from iotbx.pdb import hierarchy

import cctbx
from mmtbx import map_tools

from select_atoms import get_occupancy_groups
from utils import b_to_u_iso
from utils import chunks
from utils import get_minimum_fofc
from utils import round_step
from utils import u_iso_to_b_fac


def write_pdb_HOH_site_cart(pdb, sites_cart):
    """ Write a PDB file containing waters at supplied cartesian sites
     
    Parameters
    -------------------
    pdb: str
        Path to input PDB that contains the crystal record for the 
        crystal symmertry
    sites_cart:
        
    """

    pdb_in = hierarchy.input(file_name=pdb)

    f = open(pdb)
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

        for i, site in enumerate(sites):
            f_out.write(
                "HETATM{:>5}  O   HOH {}{:>4}{:>12.3f}{:>8.3f}{:>8.3f}  1.00 10.00           O\n".format(
                    i + 1, chain, i + 1, site[0], site[1], site[2]
                )
            )

    f_out.close()


def process_validation_csvs(start_occ, end_occ, step, set_b, out_dir, params):
    min_fofcs = []
    min_occs = []
    min_b_facs = []
    fofcs = []
    occs = []
    b_facs = []

    for lig_occupancy in np.arange(start_occ, end_occ + (step / 5), step):
        # TODO Replace CSV naming #59

        csv_name = params.exhaustive.output.csv_prefix + "_occ_{}_b_{}.csv".format(
            str(lig_occupancy).replace(".", "_"), str(set_b).replace(".", "_")
        )

        csv_path = os.path.join(out_dir, csv_name)
        min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_path)
        fofc = get_fofc_from_csv(
            csv_path, lig_occupancy, round_step(b_to_u_iso(set_b)), step
        )
        fofcs.append(fofc)
        occs.append(lig_occupancy)
        b_facs.append(set_b)
        min_b_facs.append(u_iso_to_b_fac(min_u_iso))
        min_fofcs.append(fo_fc_at_min)
        min_occs.append(min_occ)

    return min_fofcs, min_occs, min_b_facs, fofcs, occs, b_facs


def set_b_fac_all_occupancy_groups(input_pdb, output_pdb, b_fac, params):
    """ Change b factor of all atoms involved in occupancy groups"""

    pdb_inp = hierarchy.input(input_pdb)
    occ_group = get_occupancy_groups(pdb=input_pdb, params=params)
    for group in occ_group[0]:
        for residue in group:
            for chain in pdb_inp.hierarchy.only_model().chains():
                if chain.id == residue.get("chain"):
                    for residue_group in chain.residue_groups():
                        if residue_group.resseq == residue.get("resseq"):
                            for atom_group in residue_group.atom_groups():
                                for atom in atom_group.atoms():
                                    atom.set_b(b_fac)
                                    print("Changed: {}".format(residue))

    with open(output_pdb, "w") as f:
        f.write(
            pdb_inp.hierarchy.as_pdb_string(
                crystal_symmetry=pdb_inp.input.crystal_symmetry()
            )
        )


def set_b_fac_all_atoms(input_pdb, output_pdb, b_fac):
    """ Change B factors of all atoms to the same provided value"""
    pdb_inp = hierarchy.input(input_pdb)
    for chain in pdb_inp.hierarchy.only_model().chains():
        for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                for atom in atom_group.atoms():
                    atom.set_b(b_fac)

    with open(output_pdb, "w") as f:
        f.write(
            pdb_inp.hierarchy.as_pdb_string(
                crystal_symmetry=pdb_inp.input.crystal_symmetry()
            )
        )


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
        end_occ = list(set(lig_atoms))[0][1] + list(set(lig_atoms))[1][1]
        return end_occ
    else:
        "Ligand occupancy is not defined"
        exit()


def get_csv_filepath(params):
    for occupancy in np.arange(
        params.validate.options.start_simul_occ,
        params.validate.options.end_simul_occ
        + params.validate.options.step_simulation / 5,
        params.validate.options.step_simulation,
    ):
        yield os.path.join(
            params.output.out_dir,
            params.exhaustive.output.csv_prefix
            + "_occ_{}_b_{}.csv".format(
                str(occupancy).replace(".", "_"),
                str(params.validate.options.set_b).replace(".", "_"),
            ),
        )


def get_random_starting_occ_from_folder_name(occupancy, out_path, dataset_prefix):
    """Pull occupancy values from folder names of random refinements"""

    folders = [
        name
        for name in os.listdir(
            os.path.join(
                out_path,
                dataset_prefix + "_refine_occ_" + str(occupancy).replace(".", "_"),
            )
        )
        if os.path.isdir(
            os.path.join(
                out_path,
                dataset_prefix + "_refine_occ_" + str(occupancy).replace(".", "_"),
                name,
            )
        )
    ]

    folders = [name for name in folders if "exhaustive" not in name]

    for folder in folders:
        yield float("0." + folder.split("_")[-1])


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

    occ_b_df = pd.DataFrame(lig_occ_b, columns=["Atom", "Occupancy", "B_factor"])

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


def get_occ_b(refinement_dir, lig_chain, pdb_name="refine.split.bound_state.pdb"):
    """ Given a folder get occupancies & B factor of ligand in chain"""

    pdbs = get_pdbs(refinement_dir, pdb_name)

    occ_b = []

    for pdb in pdbs:

        dataset, _ = os.path.split(pdb)
        dataset = dataset.split("/")[-1]

        occ_b_df = read_ligand_occupancy_b(pdb, lig_chain)
        mean_ligand_b_factor = occ_b_df["B_factor"].mean()
        std_ligand_b_factor = occ_b_df["B_factor"].std()

        if occ_b_df.apply(lambda x: x.nunique())[1] == 1:
            lig_occ = occ_b_df.loc("Occupancy")[0][1]
            occ_b.append([dataset, lig_occ, mean_ligand_b_factor, std_ligand_b_factor])
        else:
            print(occ_b_df)
            print(
                "Occupancy varies across ligand, " "histogram not currently generated"
            )

    print(occ_b)
    occ_df = pd.DataFrame(
        data=occ_b, columns=["dataset", "occupancy", "mean_b_fac", "std_b_fac"]
    )

    return occ_df


def get_fofc_from_csv(csv_name, occupancy, u_iso, step=0.05):
    occupancy = round_step(occupancy, base=step)
    u_iso = round_step(u_iso, base=step)

    # TODO Remove this dual .csv by cleaning up csv name: issue 59

    if csv_name.endswith(".csv"):
        data = np.genfromtxt(csv_name, delimiter=",", skip_header=0)
    else:
        data = np.genfromtxt("{}.csv".format(csv_name), delimiter=",", skip_header=0)

    e = 0.0001
    data_line = data[
        ((occupancy - e) < data[:, 0])
        & ((occupancy + e) > data[:, 0])
        & ((u_iso - e) < data[:, 2])
        & ((u_iso + e) > data[:, 2])
    ]

    fo_fc = data_line[0][3]

    return fo_fc


def datasets_from_compound(protein_prefix, compound_folder):
    """ Loop over datasets with prefix in folder."""

    for dataset in filter(
        lambda x: x.startswith(protein_prefix)
        and os.path.isdir(os.path.join(compound_folder, x)),
        os.listdir(compound_folder),
    ):
        yield dataset


def collate_edstats_scores(protein_prefix, compound_folder):
    """ Collate edstats scores into DataFrame. Write csv. Return df"""

    dfs = []
    for dataset in datasets_from_compound(protein_prefix, compound_folder):

        edstats_csv = os.path.join(
            compound_folder, dataset, "Plots", "residue_scores.csv"
        )

        if os.path.exists(edstats_csv):
            edstats_df = pd.read_csv(edstats_csv)
            edstats_df["Dataset"] = dataset
            dfs.append(edstats_df)
        else:
            print("No edstats scores found for {}".format(dataset))
    try:
        compound_edstats = pd.concat(dfs, ignore_index=True)
    except ValueError:
        print("Edstats not able to concatanate")
        return None

    print(compound_edstats)

    compound_edstats.to_csv(
        file_name=os.path.join(compound_folder, "collated_edstats"), index=False
    )

    return compound_edstats


def remove_residues(input_pdb, output_pdb, residues_remove):
    pdb_in = hierarchy.input(file_name=input_pdb)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()

    selection_string_list = []
    for residue_remove in residues_remove:
        selection_string = "(resid {} and chain {} and altid {})".format(
            residue_remove[0], residue_remove[1], residue_remove[2]
        )
        selection_string_list.append(selection_string)

    selection_string = " or ".join(selection_string_list)
    not_selection_string = "not ({})".format(selection_string)

    acceptor_hierarchy = pdb_in.construct_hierarchy()

    print(not_selection_string)

    remove_atoms_sel = sel_cache.selection(not_selection_string)
    removed_hier = acceptor_hierarchy.select(remove_atoms_sel)

    f = open(os.path.join(output_pdb), "w+")

    f.write(
        removed_hier.as_pdb_string(crystal_symmetry=pdb_in.input.crystal_symmetry())
    )

    f.close()


def is_almost_equal(x, y, epsilon=1 * 10 ** (-8)):
    """Return True if two values are close in numeric value
        By default close is withing 1*10^-8 of each other
        i.e. 0.00000001
    """
    return abs(x - y) <= epsilon


def print_hier_atoms(hier):
    """
    Basic printing of an iotbx.pdb.hierarchy

    :param hier:
    :return:
    """

    for model in hier.models():
        print("Model: {}".format(model.id))

        for chain in model.chains():
            print("Chain: {}".format(chain.id))

            for residue_group in chain.residue_groups():
                print("Residue: {}".format(residue_group.resseq))

                for atom_group in residue_group.atom_groups():
                    print("Altloc: {}".format(atom_group.altloc))

                    for atom in atom_group.atoms():
                        print("Atom Name: {}".format(atom.name))

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


def compute_maps(fmodel, crystal_gridding, map_type):
    """Compute electron density maps for a given model.

    Given a model through fmodel, a map type:
    "mFo-DFc"
    "2mFo-DFc"
    Calculate a map.

    Volume scaling is applied to the map

    Return the fft map, real map, and map coefficents.

    Parameters
    ----------
    fmodel: mmtbx.f_model.f_model.manager
        cctbx object handling the model

    crystal_gridding: cctbx.maptbx.crystal_gridding
        cctbx object handling the grid on which the maps are defined

    map_type: str
        "mFo-DFc" or "2mFo-DFc" defining the map type

    Returns
    -------
    fft_map: cctbx.miller.fft_map
        Container for an FFT from reciprocal space (complex double) into real space.

    fft_map.real_map_unpadded(): scitbx_array_family_flex_ext.double
        Real component of the FFT'd map,
        removing any padding required for the FFT grid.

    map_coefficents:cctbx.miller.array object
        coeffiecients

    """

    map_coefficients = map_tools.electron_density_map(fmodel=fmodel).map_coefficients(
        map_type=map_type, isotropize=True, fill_missing=False
    )

    fft_map = cctbx.miller.fft_map(
        crystal_gridding=crystal_gridding, fourier_coefficients=map_coefficients
    )

    fft_map.apply_volume_scaling()

    return fft_map, fft_map.real_map_unpadded(), map_coefficients


def get_mean_fofc_over_cart_sites(sites_cart, fofc_map, inputs):
    """Get mean of |Fo-Fc| over a cartesian point list.
{
    Parameters
    -----------
    sites_cart: scitbx_array_family_flex_ext.vec3_double
        Cartesian sites over which to calculate the mean of |Fo-Fc|

    fofc_map: scitbx_array_family_flex_ext.double
        Real component of the FFT'd F_obs - F_calc map,
        removing any padding required for the FFT grid.

    inputs:mmtbx.utils.process_command_line_args
        holds arguments to be used for the xtal model

    Returns
    -------
    mean_abs_fofc_value: float
        Mean value of the |Fo-Fc| map over the supplied cartesian site{s

    """

    # Set a default value of parameter to sum over
    sum_abs_fofc_value = 0

    # Loop over all cartesian points
    for site_cart in list(sites_cart):
        # Get the fractional site from the cartesian coordinate
        site_frac = inputs.crystal_symmetry.unit_cell().fractionalize(site_cart)

        # Use interpolation to get the difference map value at the site
        fofc_value = fofc_map.eight_point_interpolation(site_frac)

        # Append value to sum over points
        sum_abs_fofc_value += abs(fofc_value)

    # Get the mean value of |Fo-Fc|
    mean_abs_fofc_value = sum_abs_fofc_value / len(list(sites_cart))

    return mean_abs_fofc_value


def iter_u_iso_occ(params):
    """Get occupancy and u_iso from minima, maxima and step size

    Parameters
    ----------
    params: libtbx.phil.scope_extract
            python object from phil file,
            edited with any additional parameters

    Returns
    -------
    u_iso_occ: list
        list of tuples of u_iso and occupancy to be iterated over
    """

    u_iso_occ = []
    for occupancy in np.arange(
        params.exhaustive.options.lower_occ,
        params.exhaustive.options.upper_occ + params.exhaustive.options.step / 5,
        params.exhaustive.options.step,
    ):

        for u_iso in np.arange(
            params.exhaustive.options.lower_u_iso,
            params.exhaustive.options.upper_u_iso + params.exhaustive.options.step / 5,
            params.exhaustive.options.step,
        ):
            u_iso_occ.append((occupancy, u_iso))

    return u_iso_occ
