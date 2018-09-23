from __future__ import division
from __future__ import print_function

import os
from cStringIO import StringIO

import mmtbx.utils
import iotbx.pdb
import numpy as np

from scitbx.array_family import flex
from cctbx import maptbx
from iotbx import reflection_file_utils
from phil import master_phil
from mmtbx.utils import data_and_flags_master_params

import cctbx.eltbx

def dump(obj):
  for attr in dir(obj):
    print("obj.%s = %r" % (attr, getattr(obj, attr)))

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

print(len(d_spacings))
print(len(s_sqrd))

for i in xrange(0,100):
    print(d_spacings[i], s_sqrd[i])

print(fmodel.f_obs_work().d_max_min())

print(xrs.scattering_type_registry())

dump(xrs.scattering_type_registry())

print("----------------------------------------------")

xrs.scattering_type_registry().show_summary()

print("----------------------------------------------")

print(xrs.scattering_dictionary_as_string())

print("----------------------------------------------")

dump(xrs)

print("----------------------------------------------")

dump(xrs.structure_factors)

print("----------------------------------------------")

print(xrs.scattering_type_registry().last_table())

print("----------------------------------------------")

print(cctbx.eltbx.xray_scattering.best_approximation("O"))
dump(cctbx.eltbx.xray_scattering.best_approximation("O"))

print("----------------------------------------------")

print(cctbx.eltbx.xray_scattering.best_approximation("O").show())

# Below is incorrect see https://onlinelibrary.wiley.com/doi/pdf/10.1107/S0021889801017824
# https://github.com/cctbx/cctbx_project/blob/master/cctbx/eltbx/xray_scattering/__init__.py
# and look at the
for stol_sq in fmodel.f_obs_work().sin_theta_over_lambda_sq().data():

    f0 = cctbx.eltbx.xray_scattering.best_approximation("O").at_stol_sq(stol_sq)

    occ = 0.9
    R_list = []
    for B in xrange(0,100):
        R = (occ * f0 - f0*np.exp(B*stol_sq))/(occ * f0)
        R_list.append(R)

    print(R_list)
#print(fmodel.f_obs_work().structure_factors_from_scatterers(xray_structure=xrs))

# print(fmodel.f_obs_work().structure_factors_from_scatterers(
#     xray_structure=xrs).f_calc().size())
