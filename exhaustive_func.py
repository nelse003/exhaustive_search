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


class OccBLoopCaller(object):
    """Class allowing unpickable objects to passed to loop.

    This class handles the calling of main loop,
    such that only the iterable (occupancy and b factor) changes.

    This handling class is used as some parameters are unpickable.
    These parameters must stay the same between iterations of the loop.
    This is following the example in libtx.easy_mp documentation.
    However teh easy_mp parallelisation is currently failing

    Attributes
    ----------
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

    """

    def __init__(self, xtal_model_data=None, **kwargs):
        """Initialise parameters for calculating mean |Fo-Fc|.

        Either get from instance of XtalModelData,
        or from kwargs
        """
        if xtal_model_data is not None:
            self.xrs = xtal_model_data.xrs
            self.sites_frac = xtal_model_data.sites_frac
            self.fmodel = xtal_model_data.fmodel
            self.crystal_gridding = xtal_model_data.crystal_gridding
            self.inputs = xtal_model_data.inputs
            self.params = xtal_model_data.params
            self.bound_states = xtal_model_data.bound_states
            self.ground_states = xtal_model_data.ground_states
            self.cart_points = xtal_model_data.cart_points
        else:
            for key, value in kwargs.items():
                setattr(self, key, value)

    def __call__(self, u_iso_occ):
        """
        Calculate mean|Fo-Fc| of ground/bound states with occupancy and B factor.

        Take a two item list, u_iso_occ that contains occupancy and B factor.
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
        ----------
        u_iso_occ: tuple
            occupancy and U-iso fo current iteration, stored in a tuple

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

        # Extract occupancy and u_iso
        bound_occupancy = u_iso_occ[0]
        ground_occupancy = 1 - u_iso_occ[0]
        u_iso = u_iso_occ[1]

        # make local copies of fmodel and xrs
        xrs_dc = self.xrs.deep_copy_scatterers()
        f_model_dc = self.fmodel.deep_copy()

        # Update occupancy and U_iso in local copy of xrs for ground states
        xrs_dc = self.update_xrs_over_states(
            xrs=xrs_dc,
            states=self.ground_states,
            occupancy=ground_occupancy,
            u_iso=u_iso,
        )

        # Update occupancy and U_iso in local copy of xrs for bound states
        xrs_dc = self.update_xrs_over_states(
            xrs=xrs_dc, states=self.bound_states, occupancy=bound_occupancy, u_iso=u_iso
        )

        # update local f_model with
        f_model_dc.update_xray_structure(xray_structure=xrs_dc, update_f_calc=True)

        # Generate difference map at supplied occupancy
        # and b factor for edited states
        fft_map, fofc_map, fofc = compute_maps(
            fmodel=f_model_dc,
            crystal_gridding=self.crystal_gridding,
            map_type="mFo-DFc",
        )

        # Calculate the mean value of map over supplied cart_points
        mean_abs_fofc_value = get_mean_fofc_over_cart_sites(
            self.cart_points, fofc_map, self.inputs
        )

        if self.params.exhaustive.options.generate_mtz:

            output_mtz = "testing_{}_{}.mtz".format(
                str(bound_occupancy).replace(".", "_"), str(u_iso).replace(".", "_")
            )

            if self.params.exhaustive.options.mtz_prefix is not None:
                output_mtz = self.params.exhaustive.options.mtz_prefix + output_mtz

            mtz_dataset = fofc.as_mtz_dataset(column_root_label="FOFCWT")
            mtz_object = mtz_dataset.mtz_object()
            mtz_object.write(file_name=output_mtz)

            if self.params.exhaustive.options.generate_map:
                mtz2map_args = [output_mtz]
                mtz2map(args=mtz2map_args)

        return [bound_occupancy, ground_occupancy, u_iso, mean_abs_fofc_value]

    def update_xrs_over_states(self, xrs, states, occupancy, u_iso):
        """ Update local copy of x ray structure with occupancy and u_iso

        Works for bound and ground states separately

        Parameters
        -----------
        xrs: cctbx.xray.structure.structure
            X-ray structure of interest:
            A class to describe and handle information related to a crystal structure.

        states: lst
            list of <scitbx_array_family_flex_ext.bool> objects,
            describes which residues are involded either
            bound or ground state

        occupancy: float
            occupancy for the supplied state

        u_iso: float
                U_iso at supplied occupancy


        """
        for state in states:
            num_altlocs = state[1]
            set_occupancy = occupancy / num_altlocs

            for i, site_frac in enumerate(self.sites_frac):
                if state[0][i]:
                    logger.debug(
                        "set_occ: {} occ: "
                        "{} num_altlocs: {} site_frac {}".format(
                            set_occupancy, occupancy, num_altlocs, site_frac
                        )
                    )
                    xrs.scatterers()[i].occupancy = set_occupancy
                    xrs.scatterers()[i].u_iso = u_iso

        return xrs
