from exhaustive import run as exhaustive
from phil import master_phil
import numpy as np

params =  master_phil.extract()
params.input.xtal_name = "NUDT7A-x1991"

params.input.pdb = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_search_data/covalent_ratios/NUDT7A-x1991/refine.pdb"
params.input.mtz = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_search_data/covalent_ratios/NUDT7A-x1991/refine.mtz"

params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/convex_buffer_tests"
params.settings.processes = 1

for buffer in np.arange(0,4,0.5):

    params.exhaustive.options.buffer = buffer
    params.exhaustive.output.csv_name = "{}_convex_hull_buffer.csv".format(
        str(params.exhaustive.options.buffer).replace('.','_'))
    exhaustive(params=params)

params.exhaustive.output.csv_name = "convex_hull_no_buffer.csv"
params.exhaustive.options.convex_hull_buffer = False
exhaustive(params=params)

params.exhaustive.output.csv_name = "no_convex_hull.csv"
params.exhaustive.options.convex_hull = False
exhaustive(params=params)

for buffer in np.arange(0,6,0.5):

    convex_hull_ignore_nearest = False
    params.exhaustive.options.buffer = buffer
    params.exhaustive.output.csv_name = "{}_convex_hull_buffer_ignore_nearest.csv".format(
        str(params.exhaustive.options.buffer).replace('.','_'))
    exhaustive(params=params)