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
from libtbx import easy_mp
from select_occupancy_groups import get_altloc_hier, from_altloc_generate_altloc_selection, \
                                    get_coincident_altlocs, get_occupancy_groups


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
    out_dir = "DCP2BA"
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
    processes = 8
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

def get_occupancy_group_grid_points(pdb, params):

    all_occupancy_group_cart_points = get_list_occupancy_group_grid_points(pdb, params)

    occupancy_group_cart_points = flex.vec3_double()
    for occupancy_group_points in all_occupancy_group_cart_points:
        occupancy_group_cart_points = occupancy_group_cart_points.concatenate(occupancy_group_points)

    return occupancy_group_cart_points

def get_list_occupancy_group_grid_points(pdb, params):

    """Return grid points for each occupancy group in a list"""

    all_occupancy_group_cart_points_list = []
    # Get grid points per residue basis, as residues could be seperated in space
    for altloc_hier in get_altloc_hier(pdb):

        print(pdb)

        occupancy_group_cart_points = flex.vec3_double()

        for chain in altloc_hier.only_model().chains():
            for residue_group in chain.residue_groups():

                # TODO Check: Using extract xyz instead of xrs (will this be sufficent?)
                sites_residue_cart = residue_group.atoms().extract_xyz()
                grid_min = flex.double([s - params.options.buffer for s in sites_residue_cart.min()])
                grid_max = flex.double([s + params.options.buffer for s in sites_residue_cart.max()])

                grid_residue = grid.Grid(grid_spacing = params.options.grid_spacing,
                              origin = tuple(grid_min),
                              approx_max = tuple(grid_max))

                occupancy_group_cart_points = occupancy_group_cart_points.concatenate(grid_residue.cart_points())

        print(len(occupancy_group_cart_points))
        all_occupancy_group_cart_points_list.append(occupancy_group_cart_points)

    return all_occupancy_group_cart_points_list

def get_mean_fofc_over_cart_sites(sites_cart, fofc_map, inputs):

    sum_abs_fofc_value = 0

    for site_cart in list(sites_cart):
        site_frac = inputs.crystal_symmetry.unit_cell().fractionalize(site_cart)
        fofc_value = fofc_map.eight_point_interpolation(site_frac)
        sum_abs_fofc_value += abs(fofc_value)

    mean_abs_fofc_value = sum_abs_fofc_value / len(list(sites_cart))

    return mean_abs_fofc_value


# TODO Turn the loop statements into generator, reduce repeating code
def calculate_mean_fofc(params, protein_hier, inputs, fmodel, crystal_gridding, pdb):

    xrs = protein_hier.extract_xray_structure(crystal_symmetry=inputs.crystal_symmetry)
    sites_frac = xrs.sites_frac()

    # TODO Loop over b factor separately for multiple ligands: Is this needed?
    # TODO implement an iteratively smaller step size based on minima
    u_iso_occ = []
    for occupancy in np.arange(params.options.lower_occ, params.options.upper_occ, params.options.step):
        for u_iso in np.arange(params.options.lower_u_iso, params.options.upper_u_iso, params.options.step):
            u_iso_occ.append((occupancy,u_iso))

    occupancy_group_cart_points = get_occupancy_group_grid_points(pdb, params)
    coincident_altlocs = get_coincident_altlocs(pdb)
    occupancy_groups = get_occupancy_groups(pdb)

    occ_b_loop = occ_b_loop_caller(xrs =xrs,
                                   sites_frac = sites_frac,
                                   fmodel = fmodel,
                                   crystal_gridding = crystal_gridding,
                                   inputs = inputs,
                                   params = params ,
                                   coincident_altlocs = coincident_altlocs,
                                   occupancy_group_cart_points = occupancy_group_cart_points,
                                   occupancy_groups= occupancy_groups,
                                   pdb = pdb)

    sum_fofc_results = easy_mp.pool_map(fixed_func = occ_b_loop, args = u_iso_occ, processes = params.options.processes)

    with open(params.options.csv_name + ".csv", 'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        writer.writerows(sum_fofc_results)
        sys.stdout.flush()

class occ_b_loop_caller(object):
    """This class handles the calling of main loop, such that only the iterable (occupancy and b factor) changes.
    Needed as some parameters are unpickable. These parameters must stay the same between iterations of the loop.
    """
    def __init__(self, xrs, sites_frac, fmodel, crystal_gridding, inputs, params,coincident_altlocs,
                                    occupancy_group_cart_points, occupancy_groups, pdb):
        self.xrs = xrs
        self.sites_frac = sites_frac
        self.fmodel = fmodel
        self.crystal_gridding = crystal_gridding
        self.inputs = inputs
        self.params = params
        self.coincident_altlocs = coincident_altlocs
        self.occupancy_groups = occupancy_groups
        self.occupancy_group_cart_points = occupancy_group_cart_points
        self.pdb = pdb

    def __call__(self, u_iso_occ):
        return calculate_fofc_occupancy_b_factor(u_iso_occ,
                                                 xrs=self.xrs,
                                                 sites_frac=self.sites_frac,
                                                 fmodel=self.fmodel,
                                                 crystal_gridding=self.crystal_gridding,
                                                 inputs=self.inputs,
                                                 params=self.params,
                                                 coincident_altlocs=self.coincident_altlocs,
                                                 occupancy_groups = self.occupancy_groups,
                                                 occupancy_group_cart_points=self.occupancy_group_cart_points,
                                                 pdb=self.pdb)

def calculate_fofc_occupancy_b_factor(iter_u_iso_occ,
                                      xrs,
                                      sites_frac,
                                      fmodel,
                                      crystal_gridding,
                                      inputs,
                                      params,
                                      coincident_altlocs,
                                      occupancy_groups,
                                      occupancy_group_cart_points,
                                      pdb):

    """ Main loop over which mean fofc is calculated, given occupancy and B factor """

    occupancy = iter_u_iso_occ[0]
    u_iso = iter_u_iso_occ[1]

    xrs_dc = xrs.deep_copy_scatterers()

    for altloc_group in coincident_altlocs:
        for altloc in altloc_group:
            occupancy_to_set = occupancy / len(altloc_group)
            altloc_selection = from_altloc_generate_altloc_selection(pdb, altloc, occupancy_groups)

            # TODO Check if altloc selection here will work on protein hier,
            # TODO if not redo selection based on generated selection string

            for i, site_frac in enumerate(sites_frac):
                if (altloc_selection[i]):
                    xrs_dc.scatterers()[i].occupancy = occupancy_to_set
                    xrs_dc.scatterers()[i].u_iso = u_iso

        # Alter occupancy between altloc groups
        occupancy = 1 - occupancy

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

    mean_abs_fofc_value = get_mean_fofc_over_cart_sites(occupancy_group_cart_points, fofc_map, inputs)

    row = [occupancy, u_iso, mean_abs_fofc_value]
    return row

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

    """ Main Function, Setup for protien model and run mean |Fo-Fc| calculation. 
    
    Currently selects ligands based on chains found in split.bound.pdb bs split.ground.pdb"""

    # TODO Check how this works with giant.run_default and pick correct method for loading function and holing as command scrript
    params = master_phil.extract()

    ####################################################
    # Read input PDB and reflection files. Parse into fmodel and hierarchies
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

    print(os.getcwd())

    #####################################################
    # Organise output directory
    output_folder = "{}".format(xtal_name)
    output_path = os.path.join(os.getcwd(),output_folder)
    #output_path_base = os.path.join(os.getcwd(),"NUDT22A")

    # print(output_path_base)
    #
    # if not os.path.exists(output_path_base):
    #      os.mkdir(output_path_base)
    #
    if not os.path.exists(output_path):
         os.mkdir(output_path)
    # # TODO Swap out the change directory to write out for each function writing out.
    os.chdir(output_folder)

    print(os.getcwd())

    ####################################################

    # Choose between looping over 1 ligand structure or 2
    pdb = args[0]

    # Run main calculation of |Fo-Fc| at grid points near ligand
    calculate_mean_fofc(params = params,
                        protein_hier = ph,
                        inputs = inputs,
                        fmodel = fmodel,
                        crystal_gridding = crystal_gridding,
                        pdb = pdb)
    os.chdir("../../")

if(__name__ == "__main__"):
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, blank_arg_prepend=blank_arg_prepend, args = sys.argv[1:])