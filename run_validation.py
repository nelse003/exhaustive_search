""" Basic definiton of how to setup parameter to run the validation script.

Similiar to the test_validation script, but quicker to edit"""

import os

from exhaustive.validation.validation import run as validate
from phil import master_phil

params =  master_phil.extract()

params = master_phil.extract()
params.input.xtal_name = "FALZA-x0085"
params.input.in_path = os.path.join(os.path.realpath(
    "./exhaustive/test/resources"), params.input.xtal_name)
params.validate.input.base_mtz = os.path.join(params.input.in_path,
                                                   "FALZA-x0085.free.mtz")
params.input.mtz = os.path.join(params.input.in_path,
                                     "FALZA-x0085.free.mtz")
params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")
params.output.out_dir = os.path.realpath("./exhaustive/test/output")
params.output.log_dir = os.path.join(params.output.out_dir, "logs")
params.validate.input.ground_state_pdb_path = os.path.join(
    params.input.in_path, "refine.output.ground-state.pdb")
params.validate.input.bound_state_pdb_path = os.path.join(
    params.input.in_path, "refine.output.bound-state.pdb")
params.validate.options.set_b = 40.0
params.exhaustive.options.column_type = "FMODEL"
params.exhaustive.options.generate_mtz = False
params.validate.options.use_qsub = False
params.validate.options.step_simulation = 0.1
params.validate.options.overwrite = False
params.exhaustive.options.step = 0.02
params.settings.processes = 24

if not os.path.exists(params.output.out_dir):
    os.mkdir(params.output.out_dir)

if not os.path.exists(params.output.log_dir):
    os.mkdir(params.output.log_dir)

# Removal of existing output files for cctbx fmodel to run
if os.path.exists(params.output.out_dir):

    for item in os.listdir(params.output.out_dir):
        if item.endswith(".mtz"):
            os.remove(os.path.join(params.output.out_dir, item))


validate(params)
