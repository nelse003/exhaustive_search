import logging
from cStringIO import StringIO

from cctbx import maptbx
import mmtbx.utils
import iotbx.pdb
import mmtbx.f_model
import mmtbx.masks
from mmtbx.utils import data_and_flags_master_params
from iotbx import reflection_file_utils

from utils.log_utils import log

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