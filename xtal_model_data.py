import logging
import numpy as np
from cStringIO import StringIO

from cctbx import maptbx
import mmtbx.utils
import iotbx.pdb
import mmtbx.f_model
import mmtbx.masks
from mmtbx.utils import data_and_flags_master_params
from iotbx import reflection_file_utils

from utils.log_utils import log
from utils.select_atoms import get_bound_ground_states

from utils.convex_hull import atom_points_from_sel_string
from utils.convex_hull import convex_hull_from_states
from utils.convex_hull import convex_hull_grid_points
from utils.convex_hull import convex_hull_per_residue

logger = logging.getLogger(__name__)


class XtalModelData(object):

    """
    Class holding crystal model and data

    Processing pdb and mtz file
    using mmtbx and cctbx to get required objects
    now encapsulated in this class for exhaustive search

    Attributes
    ----------
    pdb: str
        pdb filepath

    mtz: str
        mtz filepath

    xrs: cctbx.xray.structure
        Xray structure object

    inputs: mmtbx.utils.process_command_line_args
        holder for input arguments

    fmodel: mmtbx.f_model.f_model.manager
        cctbx object handling the model

    crystal_gridding: cctbx.maptbx.crystal_gridding
        cctbx object handling the grid on which the maps are defined

    params: libtbx.phil.scope_extract
            python object from phil file,
            edited with any additional parameters

    ground_states: list
        list containing the atoms in the ground state.
        list is shaped:
        [[selection, number of altlocs],[selection, number of altlocs]...]
        where selection is a iotbx.pdb selection object of type
        scitbx_array_family_flex_ext.bool

    bound_states: list
        list containing the atoms in the ground state.
        list is shaped:
        [[selection, number of altlocs],[selection, number of altlocs]...]
        where selection is a iotbx.pdb selection object of type
        scitbx_array_family_flex_ext.bool

    _f_obs: cctbx.miller.array
        observed data

    _r_free_flags: cctbx.miller.array
        r free flags

    Methods
    -------
    _get_f_obs_r_free()
        get r_free flags and f_observed

    _get_xrs()
        get xray structure object

    _get_fmodel()
        get model

    _determine_states()
        get ground and bound states
    """

    @log(in_msg="Creating Xtal Model Data: Processing input PDB and reflection files." \
                 "Parse into xray structure, fmodel and hierarchies",
         out_msg="Created Xtal Model Data")
    def __init__(self, params):
        """
        Create XtalModelData class to hold data

        Parameters
        ----------
        params: libtbx.phil.scope_extract
            python object from phil file,
            edited with any additional parameters
        """

        self.params = params
        self.pdb = self.params.input.pdb
        self.mtz = self.params.input.mtz

        self.inputs = mmtbx.utils.process_command_line_args(args=[self.params.input.pdb,
                                                                  self.params.input.mtz])

        self.xrs = self._get_xrs()

        self._f_obs, self._r_free_flags = self._get_f_obs_r_free()

        self.crystal_gridding = self._f_obs.crystal_gridding(d_min=self._f_obs.d_min(),
                                                             symmetry_flags=maptbx.use_space_group_symmetry,
                                                             resolution_factor=1. / 4)

        self.fmodel = self._get_fmodel()

        self.ground_states, self.bound_states = self._determine_states()

    @log(in_msg="Getting Rfree flags and Fobs",
         out_msg="Succeeded getting Rfree flags and Fobs")
    def _get_f_obs_r_free(self):
        """ Get f_obs and r_free_flags from pdb and mtz via self.inputs

        Returns
        -------
        f_obs: cctbx.miller.array
            observed data

        r_free_flags: cctbx.miller.array
            r free flags
        """

        rfs = reflection_file_utils.reflection_file_server(
            crystal_symmetry=self.inputs.crystal_symmetry,
            force_symmetry=True,
            reflection_files=self.inputs.reflection_files,
            err=StringIO())

        # Parameter object descrinign the may to process data
        # From mmtbx.utils
        data_flags_params = data_and_flags_master_params().extract()
        data_flags_params.labels = self.params.exhaustive.options.column_type

        determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
            reflection_file_server=rfs,
            parameters=data_flags_params,
            keep_going=True,
            log=StringIO())

        # Extract the need parts of determined_data_and_flags
        f_obs = determined_data_and_flags.f_obs
        r_free_flags = determined_data_and_flags.r_free_flags

        return f_obs, r_free_flags


    @log(in_msg="Getting xtal structure",
         out_msg="Suceeded in getting xtal structure")
    def _get_xrs(self):
        """ Get X-ray Structure Object

        Returns
        --------
        xrs: cctbx.xray.structure
            Xray structure object
        """

        # Parse input pdb
        pdb_inp = iotbx.pdb.input(file_name=self.inputs.pdb_file_names[0])

        # Create iotbx hierarchy
        ph = pdb_inp.construct_hierarchy()

        # Convert to xray structure object
        xrs = ph.extract_xray_structure(crystal_symmetry=self.inputs.crystal_symmetry)

        return xrs


    @log(in_msg="Getting model as fmodel object",
         out_msg="Suceeded in model as fmodel object")
    def _get_fmodel(self):

        """Get fmodel object

        Returns
        -------
        fmodel: mmtbx.f_model.f_model.manager
            cctbx object handling the model
        """

        mask_params = mmtbx.masks.mask_master_params.extract()
        mask_params.ignore_hydrogens = False
        mask_params.ignore_zero_occupancy_atoms = False

        fmodel = mmtbx.f_model.manager(f_obs=self._f_obs,
                                       r_free_flags=self._r_free_flags,
                                       mask_params=mask_params,
                                       xray_structure=self.xrs)
        fmodel.update_all_scales()

        logger.info("r_work: {0} r_free: {1}".format(fmodel.r_work(), fmodel.r_free()))

        return fmodel


    def _determine_states(self):
        """Determine the ground and bound states from pdb

        Parameters
        ----------

        Returns
        -------
        ground_states:
        bound_states:
        """
        try:
            bound_states, \
            ground_states = get_bound_ground_states(self.pdb, self.params)

        except UnboundLocalError:
            logger.info("Insufficient state information for pdb file %s", pdb)
            raise

        return bound_states, ground_states


    def iter_u_iso_occ(self):
        """Get occupancy and u_iso from minima, maxima and step size

        Parameters
        ----------

        Returns
        -------
        u_iso_occ: list
            list of tuples of u_iso and occupancy to be iterated over
        """

        u_iso_occ = []
        for occupancy in np.arange(self.params.exhaustive.options.lower_occ,
                                   self.params.exhaustive.options.upper_occ
                                   + self.params.exhaustive.options.step / 5,
                                   self.params.exhaustive.options.step):

            for u_iso in np.arange(self.params.exhaustive.options.lower_u_iso,
                                   self.params.exhaustive.options.upper_u_iso
                                   + self.params.exhaustive.options.step / 5,
                                   self.params.exhaustive.options.step):
                u_iso_occ.append((occupancy, u_iso))

        return u_iso_occ


    def calculate_mean_fofc(self, pdb):
        """Generate csv of occupancy and B factor for bound and ground states.

        Wrapper to prepare for main loop.
        Outputs a csv with ground_occupancy,
        bound_occupancy, u_iso and mean(|Fo-Fc|).

        Returns
        --------
        None

        # TODO Loop over b factor separately for multiple ligands #33
        """
        sites_frac = self.xrs.sites_frac()
        u_iso_occ = self.iter_u_iso_occ()


        logger.debug("U_ISO_OCC {}".format(u_iso_occ))



        if self.params.exhaustive.options.per_residue:

            cart_points = convex_hull_per_residue(self.pdb,
                                                  bound_states,
                                                  ground_states,
                                                  self.params)

        elif self.params.exhaustive.options.convex_hull:

            cart_points = convex_hull_from_states(self.pdb,
                                                  bound_states,
                                                  ground_states,
                                                  self.params)

        elif self.params.exhaustive.options.ligand_atom_points:

            cart_points = atom_points_from_sel_string(self.pdb,
                                                      selection_string=
                                                      self.params.exhaustive.options.atom_points_sel_string)

        elif self.params.exhaustive.options.ligand_grid_points:

            atom_points = atom_points_from_sel_string(self.pdb,
                                                      selection_string= \
                                                          self.params.exhaustive.options.atom_points_sel_string)

            cart_points = convex_hull_grid_points(atom_points,
                                                  self.params)

        else:
            cart_points = get_occupancy_group_grid_points(self.pdb,
                                                          bound_states,
                                                          ground_states,
                                                          self.params)

        logger.debug(cart_points)

        # write_pdb_HOH_site_cart(pdb=self.params.input.pdb, sites_cart=cart_points)

        logger.info("Looping over occupancy, u_iso with occupancy "
                    "betweeen {} and {} in steps of {} and u_iso "
                    "between {} and {} in steps of {}: "
                    "number of iteratrions.".format(
            self.params.exhaustive.options.lower_occ,
            self.params.exhaustive.options.upper_occ,
            self.params.exhaustive.options.step,
            self.params.exhaustive.options.lower_u_iso,
            self.params.exhaustive.options.upper_u_iso,
            self.params.exhaustive.options.step,
            len(u_iso_occ)))

        occ_b_loop = OccBLoopCaller(xrs=self.xrs,
                                    sites_frac=sites_frac,
                                    fmodel=self.fmodel,
                                    crystal_gridding=self.crystal_gridding,
                                    inputs=self.inputs,
                                    params=self.params,
                                    bound_states=bound_states,
                                    ground_states=ground_states,
                                    cart_points=cart_points)

        # TODO Investigate parallelisation

        # if self.params.settings.processes > 1:

        # For covalent ratios this wasn't working at all.
        # The map method seems fast enough even at 0.01

        sum_fofc_results = easy_mp.pool_map(fixed_func=occ_b_loop, args=u_iso_occ,
                                            processes=self.params.settings.processes)
        #sum_fofc_results = map(occ_b_loop, u_iso_occ)

        logger.info("Loop finished.\n"
                    "Writing bound occupancy, ground_occupancy, u_iso, "
                    "mean |Fo-Fc| to CSV: {}".format(
            self.params.exhaustive.output.csv_name))

        with open(self.params.exhaustive.output.csv_name, 'w') as f1:

            writer = csv.writer(f1, delimiter=',', lineterminator='\n')
            writer.writerows(sum_fofc_results)
            sys.stdout.flush()


