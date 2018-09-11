"""Exhaustive search.

Take in pdb and mtz, or csv of pdb and mtz locations,
run a exhaustive search of B factor and occupancy to
determine unique minima of mean(|Fo-Fc|).
"""

#  Imports
from __future__ import division
from __future__ import print_function

import csv
import datetime
import logging
import os
import sys
from cStringIO import StringIO

import cctbx.miller
import iotbx.ccp4_map
import iotbx.pdb
import mmtbx.f_model
import mmtbx.masks
import mmtbx.utils
import numpy as np
from cctbx import maptbx
from iotbx import reflection_file_utils
from libtbx import easy_mp
from mmtbx import map_tools
from mmtbx.command_line.mtz2map import run as mtz2map
from mmtbx.utils import data_and_flags_master_params

from utils.convex_hull import convex_hull_from_states, \
    atom_points_from_sel_string, convex_hull_grid_points, \
    convex_hull_per_residue

from utils.utils import is_almost_equal
from utils.select_atoms import process_refined_pdb_bound_ground_states, get_occupancy_group_grid_points
from phil import master_phil

##############################################################
PROGRAM = 'Exhaustive Search'
DESCRIPTION = """
    Take in pdb and mtz, or csv of pdb and mtz locations,
    run a exhaustive search of B factor and occupancy to
    determine unique minima.
"""
blank_arg_prepend = {'.pdb': 'pdb=', '.mtz': 'mtz=', '.csv': 'csv='}
##############################################################


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
                        pdb, logging):

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

    logging.debug("U_ISO_OCC")
    logging.debug(u_iso_occ)

    try:
        bound_states,\
        ground_states = process_refined_pdb_bound_ground_states(pdb, params)

    except UnboundLocalError:
        logging.info("Insufficient state information for pdb file %s", pdb)
        raise

    if params.exhaustive.options.per_residue:

        cart_points = convex_hull_per_residue(pdb,
                                              bound_states,
                                              ground_states,
                                              params)

    elif params.exhaustive.options.convex_hull:

        cart_points = convex_hull_from_states(pdb,
                                              bound_states,
                                              ground_states,
                                              params)

    elif params.exhaustive.options.ligand_atom_points:

        cart_points = atom_points_from_sel_string(pdb,
                                                  selection_string =
                                                  params.exhaustive.options.atom_points_sel_string)


    elif params.exhaustive.options.ligand_grid_points:

        atom_points = atom_points_from_sel_string(pdb,
                                                  selection_string =
                                                  params.exhaustive.options.atom_points_sel_string)

        cart_points = convex_hull_grid_points(atom_points,
                                              params)

    else:
        cart_points = get_occupancy_group_grid_points(pdb,
                                                      bound_states,
                                                      ground_states,
                                                      params)


    logging.debug(cart_points)

    logging.info("Looping over occupancy, u_iso with occupancy "
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

    logging.info("Loop finished.\n"
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

    logging.debug("Number of fractional sites (atoms), {}".format(sites_frac.size()))

    bound_count_true = 0
    for bound_state in bound_states:
        num_altlocs = bound_state[1]
        set_bound_occupancy = bound_occupancy / num_altlocs

        for i, site_frac in enumerate(sites_frac):
            if (bound_state[0][i]):
                bound_count_true += 1
                logging.debug("set_bound_occ: {} bound_occ: {} num_altlocs: {} site_frac {}".format(
                    set_bound_occupancy,
                    bound_occupancy,
                    num_altlocs,
                    site_frac))
                xrs_dc.scatterers()[i].occupancy = set_bound_occupancy
                xrs_dc.scatterers()[i].u_iso = u_iso

    for ground_state in ground_states:

        num_altlocs = ground_state[1]
        set_ground_occupancy = ground_occupancy / num_altlocs

        for i, site_frac in enumerate(sites_frac):
            if (ground_state[0][i]):
                logging.debug("set_bound_occ: {} bound_occ: {} num_altlocs: {} site_frac {}".format(
                    set_bound_occupancy,
                    bound_occupancy,
                    num_altlocs,
                    site_frac))
                xrs_dc.scatterers()[i].occupancy = set_ground_occupancy
                xrs_dc.scatterers()[i].u_iso = u_iso

    fmodel.update_xray_structure(
        xray_structure=xrs_dc,
        update_f_calc=True)
    fft_map, fofc_map, fofc = compute_maps(
        fmodel=fmodel,
        crystal_gridding=crystal_gridding,
        map_type="mFo-DFc")

    # For debugging the pdb file, for correct number of altlocs at correct occupancies

    if is_almost_equal(bound_occupancy,0.25) \
            and is_almost_equal(u_iso,0.5):
        logging.debug("Bound occ {}, u_iso {}".format(bound_occupancy, u_iso))
        logging.debug(str(xrs_dc.as_pdb_file()))

    if is_almost_equal(bound_occupancy,0.25) \
            and is_almost_equal(u_iso,0.7):
        logging.debug("Bound occ {}, u_iso {}".format(bound_occupancy, u_iso))
        logging.debug(str(xrs_dc.as_pdb_file()))

    if is_almost_equal(bound_occupancy,0.5) \
            and is_almost_equal(u_iso,0.5):
        logging.debug("Bound occ {}, u_iso {}".format(bound_occupancy, u_iso))
        logging.debug(str(xrs_dc.as_pdb_file()))

    if is_almost_equal(bound_occupancy,0.95) \
            and is_almost_equal(u_iso,0.5):
        logging.debug("Bound occ {}, u_iso {}".format(bound_occupancy, u_iso))
        logging.debug(str(xrs_dc.as_pdb_file()))

    if params.exhaustive.options.generate_mtz:

        output_mtz = "testing_{}_{}.mtz".format(
            str(bound_occupancy).replace(".", "_"),
            str(u_iso).replace(".", "_"))

        if params.exhaustive.options.mtz_prefix is not None:
            output_mtz =  params.exhaustive.options.mtz_prefix + output_mtz

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

    if not os.path.exists(os.path.join(params.output.out_dir,params.output.log_dir)):
        os.mkdir(os.path.join(params.output.out_dir,params.output.log_dir))


    log_time = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M.log")
    log_path = os.path.join(params.output.out_dir,
                            params.output.log_dir,
                            params.exhaustive.output.log_name + log_time)
    logging.basicConfig(filename=log_path, level=logging.DEBUG)
    logging.info("Running Exhaustive Search \n\n")

    modified_phil = master_phil.format(python_object=params)
    logging.info("Current Parameters")
    logging.info(master_phil.format(python_object=params).as_str())
    logging.info("Parameters Different from default")
    logging.info(master_phil.fetch_diff(source=modified_phil).as_str())

    #logging = start_exhaustive_logging(params)

    args = [params.input.pdb, params.input.mtz]

    header = " ############################################# "
    logging.info("\n {} \n #".format(header) + params.input.xtal_name
                + ": running exhaustive search \n {}".format(header))

    logging.info("Processing input PDB and reflection files. "
                "Parse into xray structure, fmodel and hierarchies")

    inputs = mmtbx.utils.process_command_line_args(args=args)
    logging.debug("Processed command line arguments using mmtbx.utils")

    rfs = reflection_file_utils.reflection_file_server(
        crystal_symmetry=inputs.crystal_symmetry,
        force_symmetry=True,
        reflection_files=inputs.reflection_files,
        err=StringIO())
    logging.debug("Processed reflection files using reflection file server")

    column_type = params.exhaustive.options.column_type
    logging.debug("Extracting a copy of data_and_flags_master_params "
                 "from mmtbx utils. Adding labels {} for mtz column type "
                 "to use".format(column_type))
    data_flags_params = data_and_flags_master_params().extract()
    data_flags_params.labels = column_type

    logging.debug("Default parameters supplied to "
                 "mmtbx.utils.determine_data_and_flags")
    logging.debug(data_and_flags_master_params().as_str())

    logging.debug("Current parameters supplied to "
                 "mmtbx.utils.determine_data_and_flags")
    logging.debug(data_and_flags_master_params().format(
        python_object=data_flags_params).as_str())

    logging.debug("Parameters different to default, "
                 "supplied to mmtbx.utils.determine_data_and_flags")
    logging.debug(data_and_flags_master_params().fetch_diff(
        source=data_and_flags_master_params().format(
            python_object=data_flags_params)).as_str())

    determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
        reflection_file_server=rfs,
        parameters=data_flags_params,
        keep_going=True,
        log=StringIO())

    logging.debug("Processed data and flags")

    pdb_inp = iotbx.pdb.input(file_name=inputs.pdb_file_names[0])

    logging.debug("Constructing hierarchy from input PDB: "
                 + inputs.pdb_file_names[0])

    ph = pdb_inp.construct_hierarchy()
    xrs = ph.extract_xray_structure(
        crystal_symmetry=inputs.crystal_symmetry)

    # TODO To log as string #68
    xrs.show_summary()

    logging.info("Extract Fobs and free-r flags")

    f_obs = determined_data_and_flags.f_obs
    r_free_flags = determined_data_and_flags.r_free_flags

    logging.info("Define map grididng")

    crystal_gridding = f_obs.crystal_gridding(
        d_min=f_obs.d_min(),
        symmetry_flags=maptbx.use_space_group_symmetry,
        resolution_factor=1./4)

    logging.info("Define fmodel")

    mask_params = mmtbx.masks.mask_master_params.extract()
    mask_params.ignore_hydrogens = False
    mask_params.ignore_zero_occupancy_atoms = False
    fmodel = mmtbx.f_model.manager(
        f_obs=f_obs,
        r_free_flags=r_free_flags,
        mask_params=mask_params,
        xray_structure=xrs)
    fmodel.update_all_scales()
    logging.info("r_work: {0} r_free: {1}".format(fmodel.r_work(),
                                                 fmodel.r_free()))
    logging.info("Organising output directory")
    os.chdir(params.output.out_dir)

    pdb = args[0]

    logging.info("Run main calculation of |Fo-Fc| at grid points near ligand")

    try:
        calculate_mean_fofc(params=params,
                            xrs=xrs,
                            inputs=inputs,
                            fmodel=fmodel,
                            crystal_gridding=crystal_gridding,
                            pdb=pdb,
                            logging=logging)
    except UnboundLocalError:
        raise

    os.chdir("../../")


if(__name__ == "__main__"):
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil,
                blank_arg_prepend=blank_arg_prepend, args=sys.argv[1:])
