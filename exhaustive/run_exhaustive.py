from exhaustive import run as exhaustive
from phil import master_phil

params =  master_phil.extract()
params.input.pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/" \
                   "initial_model/NUDT7A-x0299/refine.pdb"
params.input.mtz = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/" \
                   "initial_model/NUDT7A-x0299/refine.mtz"
params.input.xtal_name = "NUDT7A-x0299"
params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                        "exhaustive_search_data/convex_hull"
params.output.log_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                        "exhaustive_search_data/convex_hull/logs"
params.settings.processes = 20

exhaustive(params=params)