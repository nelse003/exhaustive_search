"""Exhaustive search.

Take in pdb and mtz, or csv of pdb and mtz locations,
run a exhaustive search of B factor and occupancy to
determine unique minima of mean(|Fo-Fc|).
"""

#  Imports
from __future__ import division
from __future__ import print_function

import csv
import os
import sys
import numpy as np
import logging
import datetime
from cStringIO import StringIO
from utils.select_atoms import process_refined_pdb_bound_ground_states
from phil import master_phil
import cctbx.miller
import giant.grid as grid
import iotbx.ccp4_map
import iotbx.pdb
import libtbx.phil
import mmtbx.f_model
import mmtbx.masks
import mmtbx.utils
from cctbx import maptbx
from iotbx import reflection_file_utils
from libtbx import easy_mp
from mmtbx import map_tools
from scitbx.array_family import flex
from scipy.spatial import ConvexHull
from scipy.spatial import distance
from mmtbx.utils import data_and_flags_master_params
from mmtbx.command_line.mtz2map import run as mtz2map
from giant.structure.select import find_nearest_atoms, protein
from utils.select_atoms import get_occupancy_groups

##############################################################
PROGRAM = 'Exhaustive Search'
DESCRIPTION = """
    Take in pdb and mtz, or csv of pdb and mtz locations,
    run a exhaustive search of B factor and occupancy to
    determine unique minima.
"""
blank_arg_prepend = {'.pdb': 'pdb=', '.mtz': 'mtz=', '.csv': 'csv='}
##############################################################

def start_exhaustive_logger(params):
    """Prepare logging.

    Logging for exhaustive search using python logging module.
    Includes a description of parameters, and parameter different to default.
    """
    log_time = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M.log")
    log_path = os.path.join(params.output.log_dir,
                            params.exhaustive.output.log_name + log_time)
    hdlr = logging.FileHandler(log_path)
    logger = logging.getLogger(__name__)
    formatter = logging.Formatter('%(asctime)s %(levelname)s \n %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)

    logger.info("Running Exhaustive Search \n\n")

    modified_phil = master_phil.format(python_object=params)
    logger.info("Current Parameters")
    logger.info(master_phil.format(python_object=params).as_str())
    logger.info("Parameters Different from default")
    logger.info(master_phil.fetch_diff(source=modified_phil).as_str())

    return logger


def compute_maps(fmodel, crystal_gridding, map_type):
    """Compute electron density maps for a given model.

    Given a model through fmodel, a map type:
    "mFo-DFc"
    "2mFo-DFc"
    Calculate a map.
    Return the fft map, real map, and map coefficents.

    :param fmodel: fmodel is a class that contains the
    :type fmodel:
    :param crystal_gridding: Seperation of crystal grid
    :type crystal_gridding: ????????????
    :param map_type: Specify the type of map to be generate i.e mFo-DFc
    :type map_type: str
    :return: fft_map.real_map_unpadded():
    :type
    """
    map_coefficients = map_tools.electron_density_map(
        fmodel=fmodel).map_coefficients(
        map_type=map_type,
        isotropize=True,
        fill_missing=False)

    fft_map = cctbx.miller.fft_map(
        crystal_gridding=crystal_gridding,
        fourier_coefficients=map_coefficients)

    fft_map.apply_volume_scaling()

    return fft_map, fft_map.real_map_unpadded(), map_coefficients


def get_occupancy_group_grid_points(pdb, bound_states, ground_states,
                                    params, logger):
    """Produce cartesian points related to occupancy groups.

    Get cartesian points that correspond to atoms involved in the
    occupancy groups (as in multi-state.restraints.params)

    :param pdb: Input PDB file
    :type path
    :param params: Working phil parameters
    :type
    :return: occupancy_group_cart_points: The cartesian points involved
    in the bound and ground states as a list
    """
    logger.info("For all bound and ground states, "
                "select cartesian grid points for each altloc/residue \n"
                "involved in occupancy groups. A buffer of {} Angstrom \n"
                "is applied to minimal and maximal grid points,"
                "with a grid seperation of {}.\n".format(
        params.exhaustive.options.buffer,
        params.exhaustive.options.grid_spacing))

    states = bound_states + ground_states

    pdb_in = iotbx.pdb.hierarchy.input(pdb)
    pdb_atoms = pdb_in.hierarchy.atoms()

    occupancy_group_cart_points = flex.vec3_double()
    for state in states:

        selection = state[0]
        selected_atoms = pdb_atoms.select(selection)
        sites_cart = selected_atoms.extract_xyz()
        grid_min = flex.double([s - params.exhaustive.options.buffer
                                for s in sites_cart.min()])
        grid_max = flex.double([s + params.exhaustive.options.buffer
                                for s in sites_cart.max()])
        grid_from_selection = grid.Grid(
            grid_spacing=params.exhaustive.options.grid_spacing,
            origin=tuple(grid_min),
            approx_max=tuple(grid_max))

        # TODO Move to logger
        print(grid_from_selection.summary())

        occupancy_group_cart_points = occupancy_group_cart_points.concatenate(
            grid_from_selection.cart_points())

    logger.info("Number of cartesian points to calculate "
                "|Fo-Fc| over: {}".format(len(occupancy_group_cart_points)))

    return occupancy_group_cart_points

def extend_convex_hull(pdb, bound_states, ground_states, params):

    states = bound_states + ground_states
    pdb_in = iotbx.pdb.hierarchy.input(pdb)
    pdb_atoms = pdb_in.hierarchy.atoms()
    atom_points = flex.vec3_double()

    all_selected_atoms =[]

    # Select atoms that are involved in the occupancy groups
    # and remove these from the distance matrix
    # This will mean that atoms selected as close by are
    # not already in that state

    occupancy_groups = get_occupancy_groups(pdb, params)

    selection_string_list = []
    for group in occupancy_groups:
        for item in group:
            for residue in item:
                selection_string = "(chain {} " \
                                   "and altid {} " \
                                   "and resid {})".format(residue['chain'],
                                                          residue['altloc'],
                                                          residue['resseq'])
                selection_string_list.append(selection_string)

    selection_string = " or ".join(selection_string_list)
    not_selection_string ="not ({})".format(selection_string)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()
    remove_atoms_sel = sel_cache.selection(not_selection_string)
    removed_hier = pdb_in.hierarchy.select(remove_atoms_sel)

    atoms_not_in_occ_group = removed_hier.atoms()

    for state in states:

        selection = state[0]
        selected_atoms = pdb_atoms.select(selection)
        all_selected_atoms.append(selected_atoms)
        sites_cart = selected_atoms.extract_xyz()
        atom_points = atom_points.concatenate(sites_cart)

    atoms_not_in_occ_group_xyz = atoms_not_in_occ_group.extract_xyz()

    # generate a convex hull
    hull = ConvexHull(atom_points)

    # Find atoms closest to the convex hull, which are not part of the convex hull
    buffer_points = []
    for vertex in hull.vertices:

        print(list(atom_points)[vertex])
        vertex_atom_point = [list(atom_points)[vertex]]
        dist_matrix = distance.cdist(vertex_atom_point, atoms_not_in_occ_group_xyz)
        min_dist = dist_matrix.min()
        index_atom = dist_matrix.argmin()

        if min_dist == 0:
            raise ValueError("Minimal distance to nearest atom is 0.0")

        nearest_atom_point = np.array(atoms_not_in_occ_group_xyz[index_atom])

        print(nearest_atom_point)
        print(type(nearest_atom_point))
        print(vertex_atom_point)
        print(type(vertex_atom_point))

        if not params.exhaustive.options.convex_hull_ignore_nearest \
            and min_dist/2 >= params.exhaustive.options.buffer:

            buffer_point = vertex_atom_point + \
                           (params.exhaustive.options.buffer/min_dist)*(nearest_atom_point - vertex_atom_point)

            print("Buffer point at buffer distance")

        elif not params.exhaustive.options.convex_hull_ignore_nearest \
            and min_dist/2 < params.exhaustive.options.buffer:

            buffer_point = vertex_atom_point + 0.5*(nearest_atom_point - vertex_atom_point)
            print("Buffer point at mid point to nearest atom")

        else:
            buffer_point = vertex_atom_point + \
                           params.exhaustive.options.buffer * \
                           (nearest_atom_point - vertex_atom_point)/min_dist
            print("Buffer point ignores nearest atom")

        buffer_points.append(buffer_point)

    buffer_points = np.concatenate(buffer_points)

    return convex_hull_grid_points(buffer_points, params)


def convex_hull_from_occupancy_group_grid_points(pdb, bound_states,
                                                 ground_states, params, logger):

    """Mean |Fo-Fc| over a convex hull defined by bound & ground states """

    states = bound_states + ground_states
    pdb_in = iotbx.pdb.hierarchy.input(pdb)
    pdb_atoms = pdb_in.hierarchy.atoms()
    atom_points = flex.vec3_double()

    all_selected_atoms =[]
    for state in states:

        selection = state[0]
        selected_atoms = pdb_atoms.select(selection)
        all_selected_atoms.append(selected_atoms)
        sites_cart = selected_atoms.extract_xyz()
        atom_points = atom_points.concatenate(sites_cart)

    return convex_hull_grid_points(atom_points,params)

def convex_hull_grid_points(atom_points,params):

    hull = ConvexHull(atom_points)

    if  isinstance(atom_points,flex.vec3_double):
        grid_min = flex.double(atom_points.min())
        grid_max = flex.double(atom_points.max())
    elif isinstance(atom_points,np.ndarray):
        grid_min = flex.double(atom_points.min(axis=0))
        grid_max = flex.double(atom_points.max(axis=0))
    else:
        raise TypeError("Convex hull grids need either numpy array, or flex.vec3_double")

    grid_from_selection = grid.Grid(
        grid_spacing=params.exhaustive.options.grid_spacing,
        origin=tuple(grid_min),
        approx_max=tuple(grid_max))

    grid_points = grid_to_np_array(grid_from_selection)

    points_in_hull = []

    for point in grid_points:
        if point_in_hull(point, hull):
            points_in_hull.append(point)

    print("Points in hull: {}".format(len(points_in_hull)))

    return flex.vec3_double(list(np.stack(points_in_hull)))

def point_in_hull(point, hull, tolerance=1e-12):
    """
    https://stackoverflow.com/questions/16750618/
    whats-an-efficient-way-to-find-if-a-
    point-lies-in-the-convex-hull-of-a-point-cl

    """

    return all(
        (np.dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)

def grid_to_np_array(grid):

    array = np.zeros((grid.cart_points().all()[0],3))
    for i, point in enumerate(grid.cart_points()):
        array[i] = point

    return array


def get_mean_fofc_over_cart_sites(sites_cart, fofc_map, inputs):
    """Get mean of |Fo-Fc| over a cartesian point list.

    :param sites_cart:
    :param fofc_map:
    :param inputs:
    :return:
    """
    sum_abs_fofc_value = 0

    for site_cart in list(sites_cart):

        site_frac = inputs.crystal_symmetry.unit_cell().\
            fractionalize(site_cart)

        fofc_value = fofc_map.eight_point_interpolation(site_frac)

        sum_abs_fofc_value += abs(fofc_value)

    mean_abs_fofc_value = sum_abs_fofc_value / len(list(sites_cart))

    return mean_abs_fofc_value


def calculate_mean_fofc(params, xrs, inputs, fmodel, crystal_gridding,
                        pdb, logger):
    """Generate csv of occupancy and B factor for bound and ground states.

    Wrapper to prepare for main loop. Outputs a csv with ground_occupancy,
    bound_occupancy, u_iso and mean(|Fo-Fc|).
    
    :param params:
    :param protein_hier:
    :param xrs:
    :param inputs:
    :param fmodel:
    :param crystal_gridding:
    :param pdb:
    :return:
    """
    sites_frac = xrs.sites_frac()

    # TODO Loop over b factor separately for multiple ligands #33
    # TODO implement an iteratively smaller step size based on minima #66

    u_iso_occ = []
    for occupancy in np.arange(params.exhaustive.options.lower_occ,
                               params.exhaustive.options.upper_occ
                               + params.exhaustive.options.step/5,
                               params.exhaustive.options.step):

        for u_iso in np.arange(params.exhaustive.options.lower_u_iso,
                               params.exhaustive.options.upper_u_iso
                               + params.exhaustive.options.step/5,
                               params.exhaustive.options.step):

            u_iso_occ.append((occupancy, u_iso))

    try:
        bound_states,\
        ground_states = process_refined_pdb_bound_ground_states(pdb, params)
    except UnboundLocalError:
        logger.info("Insufficient state information for pdb file %s", pdb)
        logger.info("Insufficient state information for pdb file %s", pdb)
        raise

    if params.testing.testing:

        extend_convex_hull(pdb, bound_states, ground_states, params)

    if params.exhaustive.options.convex_hull:

        cart_points = convex_hull_from_occupancy_group_grid_points(pdb,
                                                                      bound_states,
                                                                      ground_states,
                                                                      params,
                                                                      logger)
    else:
        cart_points = get_occupancy_group_grid_points(pdb,
                                                      bound_states,
                                                      ground_states,
                                                      params,
                                                      logger)


    logger.debug(cart_points)

    logger.info("Looping over occupancy, u_iso with occupancy "
                "betweeen {} and {} in steps of {} and u_iso "
                "between {} and {} in steps of {}.".format(
                 params.exhaustive.options.lower_occ,
                 params.exhaustive.options.upper_occ,
                 params.exhaustive.options.step,
                 params.exhaustive.options.lower_u_iso,
                 params.exhaustive.options.upper_u_iso,
                 params.exhaustive.options.step))

    occ_b_loop = occ_b_loop_caller(xrs=xrs,
                                   sites_frac=sites_frac,
                                   fmodel=fmodel,
                                   crystal_gridding=crystal_gridding,
                                   inputs=inputs,
                                   params=params,
                                   bound_states=bound_states,
                                   ground_states=ground_states,
                                   cart_points=cart_points)

    print("Pre loop")
    print(len(u_iso_occ))

    if params.settings.processes > 1:

        # For covalent ratios this wasn't working at all. The map method seems fast enough even at 0.01
        sum_fofc_results = easy_mp.pool_map(fixed_func=occ_b_loop, args=u_iso_occ,
                                            processes=params.settings.processes)
    else:
        sum_fofc_results = map(occ_b_loop,u_iso_occ)

    logger.info("Loop finished.\n"
                "Writing bound occupancy, ground_occupancy, u_iso, "
                "mean |Fo-Fc| to CSV: {}".format(
        params.exhaustive.output.csv_name))

    with open(os.path.join(params.output.out_dir,
                           params.exhaustive.output.csv_name), 'w') as f1:

        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        writer.writerows(sum_fofc_results)
        sys.stdout.flush()


class occ_b_loop_caller(object):
    """Class allowing unpickable objects to passed to loop.

    This class handles the calling of main loop,
    such that only the iterable (occupancy and b factor) changes.
    This handling class is needed as some parameters are unpickable.
    These parameters must stay the same between iterations of the loop.
    This is following the example in libtx.easy_mp documentation.
    """

    def __init__(self, xrs, sites_frac, fmodel, crystal_gridding,
                 inputs, params, bound_states, ground_states,
                 cart_points):
        """Provide all fixed parameters for calculating mean |Fo-Fc|."""
        self.xrs = xrs
        self.sites_frac = sites_frac
        self.fmodel = fmodel
        self.crystal_gridding = crystal_gridding
        self.inputs = inputs
        self.params = params
        self.bound_states = bound_states
        self.ground_states = ground_states
        self.cart_points = cart_points

    def __call__(self, u_iso_occ):
        """Calculate mean|Fo-Fc| of bound/ground states."""
        return calculate_fofc_occupancy_b_factor(
            u_iso_occ,
            xrs=self.xrs,
            sites_frac=self.sites_frac,
            fmodel=self.fmodel,
            crystal_gridding=self.crystal_gridding,
            inputs=self.inputs,
            params=self.params,
            bound_states=self.bound_states,
            ground_states=self.ground_states,
            cart_points=self.cart_points)


def calculate_fofc_occupancy_b_factor(iter_u_iso_occ,
                                      xrs,
                                      sites_frac,
                                      fmodel,
                                      crystal_gridding,
                                      inputs,
                                      params,
                                      bound_states,
                                      ground_states,
                                      cart_points):
    """Calculate mean|Fo-Fc| of ground/bound states with occupancy and B factor.

    Take a two item list, iter_u_iso that contains occupancy and B factor.
    This is used to determine mean(|Fo-Fc|) over a series of grid points
    defined by the supplied bound and ground states.

    The occupancy and B factor of scatterers (atoms) that correspond
    to the atoms in bound or ground states are set on a copy of the
    supplied xray structure (xrs_dc).

    The fmodel is updated, and used to generate |mFo-DFc| maps.
    If params.exhaustive.options.generate_mtz generate a mtz file.
    If params.exhaustive.options.generate_map generate a ccp4 map file.

    Return the occupancy of ground and bound states, u_iso
    (converatable to isotropic B factor), and the mean value of |mFo-Fc|
    over the ground and bound states.

    The B factor is set across all atoms that are varied,
    i.e. those belong to ground/ bound states which come from occupancy groups.

    :param iter_u_iso_occ:
    :param xrs:
    :param sites_frac:
    :param fmodel:
    :param crystal_gridding:
    :param inputs:
    :param params:
    :param bound_states:
    :param ground_states:
    :param occupancy_group_cart_points:
    :return:
    """
    bound_occupancy = iter_u_iso_occ[0]
    ground_occupancy = 1 - iter_u_iso_occ[0]
    u_iso = iter_u_iso_occ[1]

    xrs_dc = xrs.deep_copy_scatterers()

    for bound_state in bound_states:
        for i, site_frac in enumerate(sites_frac):
            num_altlocs = bound_state[1]
            set_bound_occupancy = bound_occupancy / num_altlocs
            if (bound_state[0][i]):
                xrs_dc.scatterers()[i].occupancy = set_bound_occupancy
                xrs_dc.scatterers()[i].u_iso = u_iso

    for ground_state in ground_states:
        for i, site_frac in enumerate(sites_frac):
            num_altlocs = ground_state[1]
            set_ground_occupancy = ground_occupancy / num_altlocs
            if (ground_state[0][i]):
                xrs_dc.scatterers()[i].occupancy = set_ground_occupancy
                xrs_dc.scatterers()[i].u_iso = u_iso

    fmodel.update_xray_structure(
        xray_structure=xrs_dc,
        update_f_calc=True)
    fft_map, fofc_map, fofc = compute_maps(
        fmodel=fmodel,
        crystal_gridding=crystal_gridding,
        map_type="mFo-DFc")

    if params.exhaustive.options.generate_mtz:

        output_mtz = "testing_{}_{}.mtz".format(
            str(bound_occupancy).replace(".", "_"),
            str(u_iso).replace(".", "_"))

        mtz_dataset = fofc.as_mtz_dataset(column_root_label="FOFCWT")
        mtz_object = mtz_dataset.mtz_object()
        mtz_object.write(file_name=output_mtz)

    if params.exhaustive.options.generate_map \
            and os.path.exists(os.path.join(params.output.out_dir,
                                            output_mtz)):

        mtz2map_args = [output_mtz]
        mtz2map(args=mtz2map_args)

    mean_abs_fofc_value = get_mean_fofc_over_cart_sites(
        cart_points, fofc_map, inputs)

    print(bound_occupancy, ground_occupancy, u_iso, mean_abs_fofc_value)

    return [bound_occupancy, ground_occupancy, u_iso, mean_abs_fofc_value]


def run(params):
    """Load protein model and run exhaustive search.

    Load in pdb and mtz file.
    Process pdb and mtz to produce:
     fmodel,
     iotbx protein hierarchy,
     xrs: xray structure
     inputs, a object used my mmtbx utils
     crystal gridding

    These are used to determine the minima of mean(|Fo-Fc|) over grid
    points defined by a search of Occupancy and B factor for atoms that
    are part of a bound ligand, or the changed protein that the bound
    atom is nearby.

    :param args:
    :param xtal_name:
    :return:
    """

    if not os.path.exists(params.output.out_dir):
        os.mkdir(params.output.out_dir)

    if not os.path.exists(params.output.log_dir):
        os.mkdir(params.output.log_dir)

    logger = start_exhaustive_logger(params)

    args = [params.input.pdb, params.input.mtz]

    header = " ############################################# "
    logger.info("\n {} \n #".format(header) + params.input.xtal_name
                + ": running exhaustive search \n {}".format(header))

    logger.info("Processing input PDB and reflection files. "
                "Parse into xray structure, fmodel and hierarchies")

    inputs = mmtbx.utils.process_command_line_args(args=args)
    logger.debug("Processed command line arguments using mmtbx.utils")

    rfs = reflection_file_utils.reflection_file_server(
        crystal_symmetry=inputs.crystal_symmetry,
        force_symmetry=True,
        reflection_files=inputs.reflection_files,
        err=StringIO())
    logger.debug("Processed reflection files using reflection file server")

    # TODO Way to select appropriate labels? #55
    column_type = "F,SIGF"
    logger.debug("Extracting a copy of data_and_flags_master_params "
                 "from mmtbx utils. Adding labels {} for mtz column type "
                 "to use".format(column_type))
    data_flags_params = data_and_flags_master_params().extract()
    data_flags_params.labels = column_type

    logger.debug("Default parameters supplied to "
                 "mmtbx.utils.determine_data_and_flags")
    logger.debug(data_and_flags_master_params().as_str())

    logger.debug("Current parameters supplied to "
                 "mmtbx.utils.determine_data_and_flags")
    logger.debug(data_and_flags_master_params().format(
        python_object=data_flags_params).as_str())

    logger.debug("Parameters different to default, "
                 "supplied to mmtbx.utils.determine_data_and_flags")
    logger.debug(data_and_flags_master_params().fetch_diff(
        source=data_and_flags_master_params().format(
            python_object=data_flags_params)).as_str())

    determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
        reflection_file_server=rfs,
        parameters=data_flags_params,
        keep_going=True,
        log=StringIO())

    logger.debug("Processed data and flags")

    pdb_inp = iotbx.pdb.input(file_name=inputs.pdb_file_names[0])

    logger.debug("Constructing hierarchy from input PDB: "
                 + inputs.pdb_file_names[0])

    ph = pdb_inp.construct_hierarchy()
    xrs = ph.extract_xray_structure(
        crystal_symmetry=inputs.crystal_symmetry)

    # TODO To log as string #68
    xrs.show_summary()
    # summary_str = ""
    # xrs.show_summary(f = summary_str)
    # print(summary_str)

    logger.info("Extract Fobs and free-r flags")

    f_obs = determined_data_and_flags.f_obs
    r_free_flags = determined_data_and_flags.r_free_flags

    logger.info("Define map grididng")

    crystal_gridding = f_obs.crystal_gridding(
        d_min=f_obs.d_min(),
        symmetry_flags=maptbx.use_space_group_symmetry,
        resolution_factor=1./4)

    logger.info("Define fmodel")

    mask_params = mmtbx.masks.mask_master_params.extract()
    mask_params.ignore_hydrogens = False
    mask_params.ignore_zero_occupancy_atoms = False
    fmodel = mmtbx.f_model.manager(
        f_obs=f_obs,
        r_free_flags=r_free_flags,
        mask_params=mask_params,
        xray_structure=xrs)
    fmodel.update_all_scales()
    logger.info("r_work: {0} r_free: {1}".format(fmodel.r_work(),
                                                 fmodel.r_free()))
    logger.info("Organising output directory")
    os.chdir(params.output.out_dir)

    pdb = args[0]

    logger.info("Run main calculation of |Fo-Fc| at grid points near ligand")

    try:
        calculate_mean_fofc(params=params,
                            xrs=xrs,
                            inputs=inputs,
                            fmodel=fmodel,
                            crystal_gridding=crystal_gridding,
                            pdb=pdb,
                            logger=logger)
    except UnboundLocalError:
        raise

    os.chdir("../../")


if(__name__ == "__main__"):
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil,
                blank_arg_prepend=blank_arg_prepend, args=sys.argv[1:])
