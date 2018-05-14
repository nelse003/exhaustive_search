from exhaustive.validation.validation import run as validate
from exhaustive.phil import master_phil
import os

params =  master_phil.extract(master_phil)

params.input.in_path = "/dls/labxchem/data/2016/lb13385-61/processing/analysis/initial_model/FALZA-x0085"
params.input.mtz = os.path.join(params.input.in_path, "FALZA-x0085.free.mtz")
params.input.pdb = os.path.join(params.input.in_path,"refine.pdb")
params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation/" \
                     "exhaustive_search_phenix_fmodel/FALZA-x0085-test-file-usage"
params.output.log_dir = os.path.join(params.output.out_dir, "logs")
params.validate.options.set_b= 40

validate(params)
