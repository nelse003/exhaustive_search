""" Basic definiton of how to setup parameter to run the validation script.

Similiar to the test_validation script, but quicker to edit"""

import os

from exhaustive.validation.validation import run as validate
from phil import master_phil
from giant.jiffies.split_conformations import master_phil as split_phil
from giant.jiffies.split_conformations import run as split_conformations

params =  master_phil.extract()

#params.input.xtal_name = "FALZA-x0085"
# params.input.in_path = os.path.join(os.path.realpath(
#     "./exhaustive/test/resources"), params.input.xtal_name)
# params.validate.input.base_mtz = os.path.join(params.input.in_path,
#                                                    "FALZA-x0085.free.mtz")
# params.input.mtz = os.path.join(params.input.in_path,
#                                      "FALZA-x0085.free.mtz")
# params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")
#params.output.out_dir = os.path.realpath("./exhaustive/test/output")

params.input.xtal_name = "NUDT7A-x6206"
params.input.in_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios_dose/NUDT7A-x6206"
params.output.out_dir = params.input.in_path
params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")
params.input.mtz = os.path.join(params.input.in_path, "refine.mtz")

params.output.log_dir = os.path.join(params.output.out_dir, "logs")
params.validate.options.set_b = 40.0
params.exhaustive.options.column_type = "FMODEL"
params.exhaustive.options.generate_mtz = True
params.validate.options.use_qsub = False
params.validate.options.step_simulation = 0.1
params.validate.options.overwrite = True
params.exhaustive.options.step = 0.05
params.settings.processes = 14

if not os.path.exists(params.output.out_dir):
    os.mkdir(params.output.out_dir)

if not os.path.exists(params.output.log_dir):
    os.mkdir(params.output.log_dir)


params.validate.input.ground_state_pdb_path = os.path.join(
    params.input.in_path, "refine.output.ground-state.pdb")
params.validate.input.bound_state_pdb_path = os.path.join(
    params.input.in_path, "refine.output.bound-state.pdb")

if not os.path.exists(params.validate.input.ground_state_pdb_path) or params.validate.options.overwrite:
    split_params = split_phil.extract()
    split_params.input.pdb = [params.input.pdb]
    split_params.output.suffix_prefix = 'output'
    split_params.options.reset_occupancies = True
    split_conformations(split_params)


# Removal of existing output files for cctbx fmodel to run
# if os.path.exists(params.output.out_dir):
#
#     for item in os.listdir(params.output.out_dir):
#         if item.endswith(".mtz"):
#             os.remove(os.path.join(params.output.out_dir, item))

validate(params)
