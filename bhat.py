from __future__ import division
from __future__ import print_function

import os
from cStringIO import StringIO

import mmtbx.utils
import iotbx.pdb

from scitbx.array_family import flex
from cctbx import maptbx
from iotbx import reflection_file_utils
from phil import master_phil
from mmtbx.utils import data_and_flags_master_params

params = master_phil.extract()
params.input.xtal_name = "FALZA-x0085"
params.input.in_path = os.path.join(os.path.realpath(
    "./exhaustive/test/resources"), params.input.xtal_name)
params.validate.input.base_mtz = os.path.join(params.input.in_path,
                                                   "FALZA-x0085.free.mtz")
params.input.mtz = os.path.join(params.input.in_path,
                                     "FALZA-x0085.free.mtz")
params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")

args = [params.input.pdb, params.input.mtz]
inputs = mmtbx.utils.process_command_line_args(args=args)

rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry=inputs.crystal_symmetry,
    force_symmetry=True,
    reflection_files=inputs.reflection_files,
    err=StringIO())

pdb_inp = iotbx.pdb.input(file_name=inputs.pdb_file_names[0])
ph = pdb_inp.construct_hierarchy()

xrs = ph.extract_xray_structure(
    crystal_symmetry=inputs.crystal_symmetry)
xrs.show_summary()

data_flags_params = data_and_flags_master_params().extract()
data_flags_params.labels = params.exhaustive.options.column_type

determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
    reflection_file_server=rfs,
    parameters=data_flags_params,
    keep_going=True,
    log=StringIO())

f_obs = determined_data_and_flags.f_obs
r_free_flags = determined_data_and_flags.r_free_flags

sites_frac = xrs.sites_frac()

f_obs.show_summary()

mask_params = mmtbx.masks.mask_master_params.extract()
mask_params.ignore_hydrogens = False
mask_params.ignore_zero_occupancy_atoms = False

fmodel = mmtbx.f_model.manager(
    f_obs=f_obs,
    r_free_flags=r_free_flags,
    mask_params=mask_params,
    xray_structure=xrs)

d_spacings = fmodel.f_obs_work().d_spacings().data()

s_sqrd = fmodel.f_obs_work().sin_theta_over_lambda_sq().data()

print(fmodel.f_obs_work().d_max_min())


#print(fmodel.f_obs_work().structure_factors_from_scatterers(xray_structure=xrs))

# print(fmodel.f_obs_work().structure_factors_from_scatterers(
#     xray_structure=xrs).f_calc().size())
