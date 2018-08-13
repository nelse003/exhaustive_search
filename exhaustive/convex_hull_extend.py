from exhaustive import run as exhaustive
from phil import master_phil

params =  master_phil.extract()
params.testing.testing = True
params.input.xtal_name = "test"

params.input.pdb = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_search_data/covalent_ratios/NUDT7A-x1906/refine.pdb"
params.input.mtz = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_search_data/covalent_ratios/NUDT7A-x1906/refine.mtz"

params.exhaustive.options.buffer = 1.0
params.settings.processes = 1
exhaustive(params=params)