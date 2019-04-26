import os
import logging

logger = logging.getLogger(__name__)

import cctbx
from mmtbx import map_tools
from mmtbx.command_line.mtz2map import run as mtz2map


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
        site_frac = inputs.crystal_symmetry.unit_cell(). \
            fractionalize(site_cart)

        # Use interpolation to get the difference map value at the site
        fofc_value = fofc_map.eight_point_interpolation(site_frac)

        # Append value to sum over points
        sum_abs_fofc_value += abs(fofc_value)

    # Get the mean value of |Fo-Fc|
    mean_abs_fofc_value = sum_abs_fofc_value / len(list(sites_cart))

    return mean_abs_fofc_value


class OccBLoopCaller(object):
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

    Parameters
    ---------
    iter_u_iso_occ: tuple
        occupancy and U-iso fo current iteration, stored in a tuple

    xrs: cctbx.xray.structure.structure
        X-ray structure of interest:
        A class to describe and handle information related to a crystal structure.

    sites_frac: scitbx_array_family_flex_ext.vec3_double
        fractional coordinates of sites

    fmodel: mmtbx.f_model.f_model.manager
        cctbx object handling the model

    crystal_gridding: cctbx.maptbx.crystal_gridding
        cctbx object handling the grid on whihc the maps are defined

    inputs: mmtbx.utils.process_command_line_args
        holds arguments to be used for the xtal model

    params: libtbx.phil.scope_extract'
        python object from phil file

    bound_states: lst
        list of <scitbx_array_family_flex_ext.bool> objects,
        describes which residues are involded in the bound state

    ground_states: lst
        list of <scitbx_array_family_flex_ext.bool> objects,
        describes which residues are involded in the bound state

    cart_points: scitbx_array_family_flex_ext.vec3_double
        cartesian coordiantes


    Returns
    -------
    list:
        bound_occupancy: float
            bound state occupancy

        ground_occupancy: float
            ground state occupancy

        u_iso: float
            U_iso at supplied occupancy

        mean_abs_fofc_value: float
            mean absolute value of difference map at
            the the supplied u_iso and occupancy,
            across the points that are around the ligand,
            and residues involved in the bound and ground states

    """

    bound_occupancy = iter_u_iso_occ[0]
    ground_occupancy = 1 - iter_u_iso_occ[0]
    u_iso = iter_u_iso_occ[1]

    xrs_dc = xrs.deep_copy_scatterers()

    logger.debug("Number of fractional sites (atoms), {}".format(sites_frac.size()))

    bound_count_true = 0
    for bound_state in bound_states:
        num_altlocs = bound_state[1]
        set_bound_occupancy = bound_occupancy / num_altlocs

        for i, site_frac in enumerate(sites_frac):
            if (bound_state[0][i]):
                bound_count_true += 1
                logger.debug("set_bound_occ: {} bound_occ: {} num_altlocs: {} site_frac {}".format(
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
                logger.debug("set_bound_occ: {} bound_occ: {} num_altlocs: {} site_frac {}".format(
                    set_bound_occupancy,
                    bound_occupancy,
                    num_altlocs,
                    site_frac))
                xrs_dc.scatterers()[i].occupancy = set_ground_occupancy
                xrs_dc.scatterers()[i].u_iso = u_iso

    # TODO Change to copy

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

        if params.exhaustive.options.mtz_prefix is not None:
            output_mtz = params.exhaustive.options.mtz_prefix + output_mtz

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

    return [bound_occupancy, ground_occupancy, u_iso, mean_abs_fofc_value]