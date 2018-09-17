import os

import mmtbx.utils

from cctbx import maptbx
from iotbx import reflection_file_utils
from phil import master_phil

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
xrs = ph.extract_xray_structure(
    crystal_symmetry=inputs.crystal_symmetry)
xrs.show_summary()