#########################################################################
#  Imports
from __future__ import print_function
from __future__ import division
import sys
from cStringIO import StringIO
import mmtbx.f_model
import mmtbx.utils
from mmtbx import map_tools
from iotbx import reflection_file_utils
import iotbx.pdb
from cctbx import maptbx
import cctbx.miller
import mmtbx.masks
import os
import csv
from scitbx.array_family import flex
import giant.grid as grid
import numpy as np

from select_hierarchies import select_chains_to_vary, chains_select_hier


# TODO change to a params style running of file
# for parameter phil file test
import libtbx.phil
##############################################################
PROGRAM = 'Exhaustive Search'
DESCRIPTION = """
    Take in pdb and mtz, or csv of pdb and mtz locations,
    run a exhaustive search of B factor and occupancy to determine unique minima.
"""
blank_arg_prepend = {'.pdb': 'pdb=','.mtz': 'mtz=','.csv': 'csv='}
##############################################################
master_phil = libtbx.phil.parse("""
input{
    pdb = None
        .type = path
    mtz = None
        .type = path
}
output{
    out_dir = "test_runs"
        .type = str
}
options{
    lower_occ = 0.0
        .type = float
    upper_occ = 1.01
        .type = float
    step = 0.05
        .type = float
    lower_u_iso = 0.2
        .type = float
    upper_u_iso = 1.21
        .type = float
    buffer = 0
        .type = float
    csv_name = 'u_iso_occupancy_vary'
        .type = str
    grid_spacing = 0.25
        .type = float
    generate_mtz = False
        .type = bool
}
""", process_includes=True)

def compute_maps(fmodel, crystal_gridding, map_type):
    map_coefficients = map_tools.electron_density_map(
        fmodel = fmodel).map_coefficients(
        map_type         = map_type,
        isotropize       = True,
        fill_missing     = False)
    fft_map = cctbx.miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = map_coefficients)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded(), map_coefficients

def pick_atoms_to_loop_over(pdb, protein_hier, inputs, fmodel, crystal_gridding, params):

    bound_pdb_path, ground_pdb_path = get_bound_ground_pdb(pdb)

    ground_bound_chains, bound_chains, ground_chains = select_chains_to_vary(bound_pdb_path, ground_pdb_path)

    if not bound_chains:
        raise ValueError("No new chains in bound structure")

    if not ground_chains:
        ("No differing atoms near bound heteroatoms in ground structure")

    print(ground_bound_chains, bound_chains, ground_chains)

    if not ground_bound_chains:
        loop_residues_altlocs_mean_fofc(params, protein_hier, inputs, fmodel, crystal_gridding, bound_chains)
    elif not ground_chains:
        loop_residues_altlocs_mean_fofc(params, protein_hier, inputs, fmodel, crystal_gridding, bound_chains)
    else:
        loop_residues_altlocs_mean_fofc_ground_bound(params, protein_hier, inputs, fmodel, crystal_gridding,
                                                     bound_chains, ground_chains, ground_bound_chains)

def get_cartesian_grid_points_near_chains(params, inputs, hier):

    # TODO try replacement with extract xyz
    xrs_chains = hier.extract_xray_structure(crystal_symmetry=inputs.crystal_symmetry)
    #print(hier.overall_counts().as_str())
    sites_cart_chains = xrs_chains.sites_cart()

    #print(params)

    # Calculate the extent of the grid
    # TODO Test different grid expansion sizes
    grid_min = flex.double([s-params.options.buffer for s in sites_cart_chains.min()])
    grid_max = flex.double([s+params.options.buffer for s in sites_cart_chains.max()])

    # TODO Remove dependence on giant.grid
    grid_near_lig = grid.Grid(grid_spacing = params.options.grid_spacing,
                              origin = tuple(grid_min),
                              approx_max = tuple(grid_max))

    return grid_near_lig.cart_points()

def get_mean_fofc_over_cart_sites(sites_cart, fofc_map, inputs):

    sum_abs_fofc_value = 0

    for site_cart in list(sites_cart):
        site_frac = inputs.crystal_symmetry.unit_cell().fractionalize(site_cart)
        fofc_value = fofc_map.eight_point_interpolation(site_frac)
        sum_abs_fofc_value += abs(fofc_value)

    mean_abs_fofc_value = sum_abs_fofc_value / len(list(sites_cart))

    return mean_abs_fofc_value

def loop_residues_altlocs_mean_fofc(params, protein_hier, inputs, fmodel, crystal_gridding, bound_chains):

    xrs = protein_hier.extract_xray_structure(crystal_symmetry=inputs.crystal_symmetry)
    sites_frac = xrs.sites_frac()

    bound_hier, sel_bound = chains_select_hier(bound_chains, protein_hier)
    sites_cart_near_lig = get_cartesian_grid_points_near_chains(params, inputs, bound_hier)

    # TODO get residue identifier?
    # TODO remove CSV name to phil file
    with open(os.path.join(params.options.csv_name, ".csv"), 'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        # currently loop over rough occupancy range for initial testing
        # TODO Loop over b factor seperately for both ligands
        for occupancy in numpy.arange(params.options.lower_occ, params.options.upper_occ, params.options.step):
            for u_iso in numpy.arange(params.options.lower_u_iso, params.options.upper_u_iso, params.options.step):
                xrs_dc = xrs.deep_copy_scatterers()

                # TODO Set occupancy using whole residue group rather than for loop
                # Change Occupancy and B factor for all atoms in selected ligand at the same time
                for i, site_frac in enumerate(sites_frac):
                    if(sel_bound[i]):
                        xrs_dc.scatterers()[i].occupancy = occupancy
                        xrs_dc.scatterers()[i].u_iso = u_iso

                fmodel.update_xray_structure(
                    xray_structure=xrs_dc,
                    update_f_calc=True)
                fofc_map, fofc = compute_maps(
                    fmodel=fmodel,
                    crystal_gridding=crystal_gridding,
                    map_type="mFo-DFc")

                mean_abs_fofc_value = get_mean_fofc_over_cart_sites(sites_cart_near_lig, fofc_map, inputs)

                row = [occupancy, u_iso, mean_abs_fofc_value]
                writer.writerow(row)
                sys.stdout.flush()



# TODO Turn the loop statements into generator, reduce repeating code
def loop_residues_altlocs_mean_fofc_ground_bound(params, protein_hier, inputs, fmodel, crystal_gridding,
                                                 bound_chains, ground_chains, bound_ground_chains):

    xrs = protein_hier.extract_xray_structure(crystal_symmetry=inputs.crystal_symmetry)
    sites_frac = xrs.sites_frac()

    bound_hier, sel_bound = chains_select_hier(bound_chains, protein_hier)
    ground_hier, sel_ground = chains_select_hier(ground_chains, protein_hier)
    bound_ground_hier, sel_bound_ground = chains_select_hier(bound_ground_chains, protein_hier)
    xrs_bound_ground = bound_ground_hier.extract_xray_structure(crystal_symmetry=inputs.crystal_symmetry)
    sites_cart_bound_ground = xrs_bound_ground.sites_cart()

    # Calculate the extent of the grid
    # TODO Test different grid expansion sizes
    grid_min = flex.double([s-params.options.buffer for s in sites_cart_bound_ground.min()])
    grid_max = flex.double([s+params.options.buffer for s in sites_cart_bound_ground.max()])

    # TODO Remove dependence on giant.grid
    grid_near_lig = grid.Grid(grid_spacing = params.options.grid_spacing,
                              origin = tuple(grid_min),
                              approx_max = tuple(grid_max))

    sites_cart_near_bound_ground = grid_near_lig.cart_points()

    # TODO implement an iteratively smaller step size based on minima

    # TODO Swap to only generate fraction sites related to selected ligand (i.e remove if statement)
    # TODO Profile the speed of the loop
    with open(params.options.csv_name + ".csv", 'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        # currently loop over rough occupancy range for initial testing
        # TODO Loop over b factor separately for both ligands: Is this needed?

        for occupancy in np.arange(params.options.lower_occ, params.options.upper_occ, params.options.step):
            for u_iso in np.arange(params.options.lower_u_iso, params.options.upper_u_iso, params.options.step):

                xrs_dc = xrs.deep_copy_scatterers()

                # TODO Set occupancy using whole residue group rather than for loop
                # Change Occupancy and B factor for all atoms in selected ligand at the same time

                # TODO Link number of copies of ligand/ bound state to divisor in the occupancy?
                for i, site_frac in enumerate(sites_frac):
                    #print(sel_ground[i])
                    if(sel_bound[i]):
                        xrs_dc.scatterers()[i].occupancy = occupancy / 2
                        xrs_dc.scatterers()[i].u_iso = u_iso
                    if(sel_ground[i]):
                        xrs_dc.scatterers()[i].occupancy = (1 - occupancy)/2
                        xrs_dc.scatterers()[i].u_iso = u_iso
                        bound_occ = xrs_dc.scatterers()[i].occupancy

                fmodel.update_xray_structure(
                    xray_structure=xrs_dc,
                    update_f_calc=True)
                fofc_map, fofc = compute_maps(
                    fmodel=fmodel,
                    crystal_gridding=crystal_gridding,
                    map_type="mFo-DFc")

                if params.options.generate_mtz:
                    mtz_dataset = fofc.as_mtz_dataset(column_root_label="FOFCWT")
                    mtz_object = mtz_dataset.mtz_object()
                    mtz_object.write(file_name="testing_{}_{}.mtz".format(occupancy, u_iso))

                mean_abs_fofc_value = get_mean_fofc_over_cart_sites(sites_cart_near_bound_ground, fofc_map, inputs)

                # TODO Link this ligand occupancy to the number of copies of the ligand in the bound state
                lig_occ = occupancy / 2

                row = [lig_occ, bound_occ, u_iso, mean_abs_fofc_value]
                writer.writerow(row)
                sys.stdout.flush()

def get_bound_ground_pdb(refinement_pdb):

    split_bound_name = os.path.basename(refinement_pdb).rstrip('.pdb') + ".split.bound-state.pdb"
    split_ground_name = os.path.basename(refinement_pdb).rstrip('.pdb') + ".split.ground-state.pdb"
    bound_pdb_path = os.path.join(os.path.dirname(refinement_pdb),split_bound_name)
    ground_pdb_path = os.path.join(os.path.dirname(refinement_pdb), split_ground_name)

    return bound_pdb_path, ground_pdb_path

def get_minimum_fofc(csv_name):
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
    #b_iso = (8 * np.pi ** 2) * u_iso ** 2

    # if three column data

    min_index = np.argmin(fo_fc)
    return occ[min_index], u_iso[min_index]

def run(args, xtal_name):

    # TODO Check how this works with giant.run_default and pick correct method for loading function and holing as command scrript
    params = master_phil.extract()

    # Read input PDB and reflection files
    inputs = mmtbx.utils.process_command_line_args(args = args)

    reflection_files = inputs.reflection_files
    rfs = reflection_file_utils.reflection_file_server(
        crystal_symmetry = inputs.crystal_symmetry,
        force_symmetry   = True,
        reflection_files = reflection_files,
        err              = StringIO())
    determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
        reflection_file_server = rfs,
        keep_going             = True,
        log                    = StringIO())
    pdb_inp = iotbx.pdb.input(file_name = inputs.pdb_file_names[0])
    ph = pdb_inp.construct_hierarchy()
    xrs = ph.extract_xray_structure(
        crystal_symmetry = inputs.crystal_symmetry)
    # Extract Fobs and free-r flags
    f_obs = determined_data_and_flags.f_obs
    r_free_flags = determined_data_and_flags.r_free_flags
    # Define map grididng
    crystal_gridding = f_obs.crystal_gridding(
        d_min             = f_obs.d_min(),
        symmetry_flags    = maptbx.use_space_group_symmetry,
        resolution_factor = 1./4)
    # Define fmodel
    mask_params = mmtbx.masks.mask_master_params.extract()
    mask_params.ignore_hydrogens=False
    mask_params.ignore_zero_occupancy_atoms=False
    fmodel = mmtbx.f_model.manager(
        f_obs          = f_obs,
        r_free_flags   = r_free_flags,
        mask_params    = mask_params,
        xray_structure = xrs)
    fmodel.update_all_scales()
    print("r_work: {0} r_free: {1}".format(fmodel.r_work(), fmodel.r_free()))


    # Generate results in an output directory

    output_folder = "{}/{}".format(params.output.out_dir, xtal_name)
    output_path = os.path.join(os.getcwd(),output_folder)
    output_path_base = os.path.join(os.getcwd(),"test_runs")

    if not os.path.exists(output_path_base):
         os.mkdir(output_path_base)

    if not os.path.exists(output_path):
        os.mkdir(output_path)
    # TODO Swap out the change directory to write out for each function writing out.
    os.chdir(output_folder)

    pdb = args[0]
    pick_atoms_to_loop_over(pdb, ph, inputs, fmodel, crystal_gridding, params)
    os.chdir("../../")

if(__name__ == "__main__"):
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, blank_arg_prepend=blank_arg_prepend, args = sys.argv[1:])