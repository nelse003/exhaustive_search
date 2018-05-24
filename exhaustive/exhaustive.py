#########################################################################
#  Imports
from __future__ import division
from __future__ import print_function

import csv
import os
import sys
from cStringIO import StringIO

from utils.select import process_refined_pdb_bound_ground_states
from phil import master_phil

import cctbx.miller
import giant.grid as grid
import iotbx.ccp4_map
import iotbx.pdb
import libtbx.phil
import mmtbx.f_model
import mmtbx.masks
import mmtbx.utils
import numpy as np
from cctbx import maptbx
from iotbx import reflection_file_utils
from libtbx import easy_mp
from mmtbx import map_tools
from scitbx.array_family import flex

from mmtbx.utils import data_and_flags_master_params

##############################################################
PROGRAM = 'Exhaustive Search'
DESCRIPTION = """
    Take in pdb and mtz, or csv of pdb and mtz locations,
    run a exhaustive search of B factor and occupancy to determine unique minima.
"""
blank_arg_prepend = {'.pdb': 'pdb=','.mtz': 'mtz=','.csv': 'csv='}
##############################################################

#########################################################################
import logging
import datetime

def start_exhaustive_logger(params):

    log_time = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M.log")
    log_path = os.path.join(params.output.log_dir, params.exhaustive.output.log_name + log_time)
    print(log_path)
    hdlr = logging.FileHandler(log_path)
    logger = logging.getLogger(__name__)
    formatter = logging.Formatter('%(asctime)s %(levelname)s \n %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)

    logger.info("Running Exhaustive Search \n\n")

    modified_phil = master_phil.format(python_object = params)
    logger.info("Current Parameters")
    logger.info(master_phil.format(python_object =params).as_str())
    logger.info("Parameters Different from default")
    logger.info(master_phil.fetch_diff(source=modified_phil).as_str())

    return logger

#########################################################################
def compute_maps(fmodel, crystal_gridding, map_type):

    """ 
    Compute electron density maps for a given model
    
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
        fmodel = fmodel).map_coefficients(
        map_type         = map_type,
        isotropize       = True,
        fill_missing     = False)
    fft_map = cctbx.miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = map_coefficients)
    fft_map.apply_volume_scaling()
    return fft_map,fft_map.real_map_unpadded(), map_coefficients


def get_occupancy_group_grid_points(pdb, bound_states, ground_states, params, logger):
    """
    Get cartesian points that correspond to atoms involved in the occupancy groups (as in multi-state.restraints.params)
    
    :param pdb: Input PDB file
    :type path
    :param params: Working phil parameters 
    :type 
    :return: occupancy_group_cart_points: The cartesian points involved in the bound and ground states as 
    """

    logger.info("For all bound and ground states, select cartesian grid points for each altloc/residue \n" \
                "involved in occupancy groups. A buffer of {} Angstrom \n".format(params.exhaustive.options.buffer) + \
                "is applied to minimal and maximal grid points, with a grid seperation of {}.\n".format(
                    params.exhaustive.options.grid_spacing))

    states = bound_states + ground_states

    pdb_in = iotbx.pdb.hierarchy.input(pdb)
    pdb_atoms = pdb_in.hierarchy.atoms()

    occupancy_group_cart_points = flex.vec3_double()
    for state in states:

        selection = state[0]
        selected_atoms = pdb_atoms.select(selection)
        sites_cart = selected_atoms.extract_xyz()
        grid_min = flex.double([s - params.exhaustive.options.buffer for s in sites_cart.min()])
        grid_max = flex.double([s + params.exhaustive.options.buffer for s in sites_cart.max()])
        grid_from_selection = grid.Grid(grid_spacing = params.exhaustive.options.grid_spacing,
                               origin = tuple(grid_min),
                               approx_max = tuple(grid_max))
        print(grid_from_selection.summary())

        occupancy_group_cart_points = occupancy_group_cart_points.concatenate(grid_from_selection.cart_points())

    logger.info("Number of cartesian points to calculate |Fo-Fc| over: {}".format(len(occupancy_group_cart_points)))
    return occupancy_group_cart_points

def get_mean_fofc_over_cart_sites(sites_cart, fofc_map, inputs):
    """
    Given cartesian sites and an |Fo-Fc| map, find the mean value of |Fo-Fc| over the cartesian points.
    
    :param sites_cart: 
    :param fofc_map: 
    :param inputs: 
    :return: 
    """

    sum_abs_fofc_value = 0

    for site_cart in list(sites_cart):
        site_frac = inputs.crystal_symmetry.unit_cell().fractionalize(site_cart)
        fofc_value = fofc_map.eight_point_interpolation(site_frac)
        sum_abs_fofc_value += abs(fofc_value)

    mean_abs_fofc_value = sum_abs_fofc_value / len(list(sites_cart))

    return mean_abs_fofc_value


def calculate_mean_fofc(params, protein_hier, xrs, inputs, fmodel, crystal_gridding, pdb,logger):
    """
    Wrapper to prepare for main loop. Outputs a csv with ground_occupancy, bound_occupancy, u_iso and mean(|Fo-Fc|). 
    
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

    # TODO Loop over b factor separately for multiple ligands: Is this needed? #33
    # TODO implement an iteratively smaller step size based on minima #66

    u_iso_occ = []
    for occupancy in np.arange(params.exhaustive.options.lower_occ,
                               params.exhaustive.options.upper_occ + params.exhaustive.options.step/5,
                               params.exhaustive.options.step):
        for u_iso in np.arange(params.exhaustive.options.lower_u_iso,
                               params.exhaustive.options.upper_u_iso + params.exhaustive.options.step/5,
                               params.exhaustive.options.step):
            u_iso_occ.append((occupancy,u_iso))

    try:
        bound_states, ground_states = process_refined_pdb_bound_ground_states(pdb, params)
    except UnboundLocalError:
        logger.info("Insufficient state information for pdb file %s",pdb)
        raise

    occupancy_group_cart_points = get_occupancy_group_grid_points(pdb, bound_states, ground_states, params, logger)

    logger.debug(occupancy_group_cart_points)

    logger.info("Looping over occupancy, u_iso with" \
                "occupancy betweeen {} and {} in steps of {}.".format(params.exhaustive.options.lower_occ,
                                                                      params.exhaustive.options.upper_occ,
                                                                      params.exhaustive.options.step) + \
                "and u_iso between {} and {} in steps of {}.".format(params.exhaustive.options.lower_u_iso,
                                                                     params.exhaustive.options.upper_u_iso,
                                                                     params.exhaustive.options.step))

    occ_b_loop = occ_b_loop_caller(xrs =xrs,
                                   sites_frac = sites_frac,
                                   fmodel = fmodel,
                                   crystal_gridding = crystal_gridding,
                                   inputs = inputs,
                                   params = params,
                                   bound_states = bound_states,
                                   ground_states = ground_states,
                                   occupancy_group_cart_points = occupancy_group_cart_points)

    sum_fofc_results = easy_mp.pool_map(fixed_func = occ_b_loop, args = u_iso_occ, processes = params.settings.processes)

    logger.info("Loop finished.\n" \
                "Writing bound occupancy, ground_occupancy, u_iso, mean |Fo-Fc| to CSV: {}".format(
        params.exhaustive.output.csv_name))

    with open(os.path.join(params.output.out_dir,params.exhaustive.output.csv_name),'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        writer.writerows(sum_fofc_results)
        sys.stdout.flush()

class occ_b_loop_caller(object):
    """This class handles the calling of main loop, such that only the iterable (occupancy and b factor) changes.
    
    This handling class is needed as some parameters are unpickable. 
    These parameters must stay the same between iterations of the loop.
    """
    def __init__(self, xrs, sites_frac, fmodel, crystal_gridding, inputs, params, bound_states, ground_states,
                 occupancy_group_cart_points):
        self.xrs = xrs
        self.sites_frac = sites_frac
        self.fmodel = fmodel
        self.crystal_gridding = crystal_gridding
        self.inputs = inputs
        self.params = params
        self.bound_states = bound_states
        self.ground_states = ground_states
        self.occupancy_group_cart_points = occupancy_group_cart_points

    def __call__(self, u_iso_occ):
        return calculate_fofc_occupancy_b_factor(u_iso_occ,
                                                 xrs=self.xrs,
                                                 sites_frac=self.sites_frac,
                                                 fmodel=self.fmodel,
                                                 crystal_gridding=self.crystal_gridding,
                                                 inputs=self.inputs,
                                                 params=self.params,
                                                 bound_states = self.bound_states,
                                                 ground_states = self.ground_states,
                                                 occupancy_group_cart_points=self.occupancy_group_cart_points)

def calculate_fofc_occupancy_b_factor(iter_u_iso_occ,
                                      xrs,
                                      sites_frac,
                                      fmodel,
                                      crystal_gridding,
                                      inputs,
                                      params,
                                      bound_states,
                                      ground_states,
                                      occupancy_group_cart_points):
    """
    Iteration of main loop over which mean(|Fo-Fc|) is calculated, given occupancy and B factor.
    
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
    ground_occupancy = 1- iter_u_iso_occ[0]
    u_iso = iter_u_iso_occ[1]

    xrs_dc = xrs.deep_copy_scatterers()

    for bound_state in bound_states:
        for i, site_frac in enumerate(sites_frac):
            num_altlocs = bound_state[1]
            set_bound_occupancy = bound_occupancy / num_altlocs
            if (bound_state[0][i]):
                xrs_dc.scatterers()[i].occupancy = set_bound_occupancy
                xrs_dc.scatterers()[i].u_iso = u_iso
                #print("gs: {}".format(xrs_dc.scatterers()[i].u_iso))

    for ground_state in ground_states:
        for i, site_frac in enumerate(sites_frac):
            num_altlocs = ground_state[1]
            set_ground_occupancy = ground_occupancy / num_altlocs
            if (ground_state[0][i]):
                xrs_dc.scatterers()[i].occupancy = set_ground_occupancy
                xrs_dc.scatterers()[i].u_iso = u_iso
                #print("bs: {}".format(xrs_dc.scatterers()[i].u_iso))

    fmodel.update_xray_structure(
        xray_structure=xrs_dc,
        update_f_calc=True)
    fft_map, fofc_map, fofc = compute_maps(
        fmodel=fmodel,
        crystal_gridding=crystal_gridding,
        map_type="mFo-DFc")

    if params.exhaustive.options.generate_mtz:
        mtz_dataset = fofc.as_mtz_dataset(column_root_label="FOFCWT")
        mtz_object = mtz_dataset.mtz_object()
        mtz_object.write(file_name="testing_{}_{}.mtz".format(str(bound_occupancy).repalce(".","_"),
                                                              str(u_iso).replace(".","_")))
    if params.exhaustive.options.generate_map:
        os.system("phenix.mtz2map testing_{}_{}.mtz".format(str(bound_occupancy).repalce(".","_"),
                                                              str(u_iso).replace(".","_")))

    # print(type(fofc_map))
    # print(type(fofc))

    mean_abs_fofc_value = get_mean_fofc_over_cart_sites(occupancy_group_cart_points, fofc_map, inputs)

    return [bound_occupancy, ground_occupancy, u_iso, mean_abs_fofc_value]

def run(params):

# TODO Clean up ability to run as command line, and source as python #67
#def run(args, xtal_name)

    """
    Main Function, Setup for protein model and run mean |Fo-Fc| calculation. 
    
    Currently selects ligands based on chains found in split.bound.pdb bs split.ground.pdb
    
    :param args: 
    :param xtal_name: 
    :return: 
    """

    ####################################################
    logger = start_exhaustive_logger(params)

    args = [params.input.pdb,params.input.mtz]

    header = " ############################################# "
    logger.info("\n {} \n #".format(header) + params.input.xtal_name+ ": running exhaustive search \n {}".format(header))

    logger.info("Processing input PDB and reflection files. Parse into xray structure, fmodel and hierarchies")

    inputs = mmtbx.utils.process_command_line_args(args = args)
    logger.debug("Processed command line arguments using mmtbx.utils")

    rfs = reflection_file_utils.reflection_file_server(
        crystal_symmetry = inputs.crystal_symmetry,
        force_symmetry   = True,
        reflection_files = inputs.reflection_files,
        err              = StringIO())
    logger.debug("Processed reflection files using reflection file server")

    # TODO Way to select appropriate labels? #55
    column_type = "F,SIGF"
    logger.debug("Extracting a copy of data_and_flags_master_params from mmtbx utils \n"
                 "Adding labels {} for mtz column type to use".format(column_type))
    data_flags_params = data_and_flags_master_params().extract()
    data_flags_params.labels=column_type
    logger.debug("Default parameters supplied to mmtbx.utils.determine_data_and_flags")
    logger.debug(data_and_flags_master_params().as_str())
    logger.debug(" Current parameters supplied to mmtbx.utils.determine_data_and_flags")
    logger.debug(data_and_flags_master_params().format(python_object =data_flags_params).as_str())
    logger.debug(" Parameters different to default, supplied to mmtbx.utils.determine_data_and_flags")
    logger.debug(data_and_flags_master_params().fetch_diff(
        source = data_and_flags_master_params().format(
            python_object =data_flags_params)).as_str())

    determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
        reflection_file_server = rfs,
        parameters = data_flags_params,
        keep_going             = True,
        log                    = StringIO())
    logger.debug("Processed data and flags")

    pdb_inp = iotbx.pdb.input(file_name = inputs.pdb_file_names[0])

    logger.debug("Constructing hierarchy from input PDB: " + inputs.pdb_file_names[0])

    ph = pdb_inp.construct_hierarchy()
    xrs = ph.extract_xray_structure(
        crystal_symmetry = inputs.crystal_symmetry)

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
        d_min             = f_obs.d_min(),
        symmetry_flags    = maptbx.use_space_group_symmetry,
        resolution_factor = 1./4)

    logger.info("Define fmodel")

    mask_params = mmtbx.masks.mask_master_params.extract()
    mask_params.ignore_hydrogens=False
    mask_params.ignore_zero_occupancy_atoms=False
    fmodel = mmtbx.f_model.manager(
        f_obs = f_obs,
        r_free_flags   = r_free_flags,
        mask_params    = mask_params,
        xray_structure = xrs)
    fmodel.update_all_scales()
    logger.info("r_work: {0} r_free: {1}".format(fmodel.r_work(), fmodel.r_free()))

    logger.info("Organising output directory")

    # output_folder = "{}".format(xtal_name)
    # output_path = os.path.join(os.getcwd(),output_folder)
    # print("OUT:{}".format(output_path))
    #output_path_base = os.path.join(os.getcwd(),"NUDT22A")

    # print(output_path_base)
    #
    # if not os.path.exists(output_path_base):
    #      os.mkdir(output_path_base)
    #
    # if not os.path.exists(output_path):
    #      os.mkdir(output_path)
    os.chdir(params.output.out_dir)

    pdb = args[0]

    logger.info("Run main calculation of |Fo-Fc| at grid points near ligand")

    try:
        calculate_mean_fofc(params = params,
                        protein_hier = ph,
                        xrs=xrs,
                        inputs = inputs,
                        fmodel = fmodel,
                        crystal_gridding = crystal_gridding,
                        pdb = pdb,
                        logger= logger)
    except UnboundLocalError:
        raise

    exit()
    os.chdir("../../")

if(__name__ == "__main__"):
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, blank_arg_prepend=blank_arg_prepend, args = sys.argv[1:])
