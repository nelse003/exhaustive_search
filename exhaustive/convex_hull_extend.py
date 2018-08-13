from exhaustive import run as exhaustive
from phil import master_phil
from plotting.plot import scatter_plot
import numpy as np
import os

scatter_plot("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/" \
"repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/" \
"FMOPL000622a_DSI_poised/NUDT22A-x0955/NUDT22A-x0955_exhaustive_search_occ_u_iso.csv",
three_dim_plot=True)

exit()

params =  master_phil.extract()
params.input.xtal_name = "NUDT7A-x1991"

params.input.pdb = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_search_data/covalent_ratios/NUDT7A-x1991/refine.pdb"
params.input.mtz = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_search_data/covalent_ratios/NUDT7A-x1991/refine.mtz"

params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/convex_buffer_tests"
params.settings.processes = 1

params.exhaustive.options.buffer = 1.35
params.exhaustive.output.csv_name = "{}_convex_hull_buffer.csv".format(
    str(params.exhaustive.options.buffer).replace('.', '_'))

#exhaustive(params=params)
scatter_plot(os.path.join(params.output.out_dir,
                          params.exhaustive.output.csv_name),
             three_dim_plot=True,
             title_text="{}_convex_hull_buffer".format(str(
                 params.exhaustive.options.buffer).replace('.', '_')))

for buffer in np.arange(0,4,0.5):

    params.exhaustive.options.buffer = buffer
    params.exhaustive.output.csv_name = "{}_convex_hull_buffer.csv".format(
        str(params.exhaustive.options.buffer).replace('.','_'))
    #exhaustive(params=params)
    scatter_plot(os.path.join(params.output.out_dir,
                          params.exhaustive.output.csv_name),
                              three_dim_plot=True,
                 title_text="{}_convex_hull_buffer".format(str(
                     params.exhaustive.options.buffer).replace('.', '_')))

scatter_plot(os.path.join(params.output.out_dir,
                          params.exhaustive.output.csv_name),
             three_dim_plot=True,
             title_text="convex_hull_no_buffer")

params.exhaustive.output.csv_name = "no_convex_hull.csv"
params.exhaustive.options.convex_hull = False
scatter_plot(os.path.join(params.output.out_dir,
                          params.exhaustive.output.csv_name),
             three_dim_plot=True,
             title_text="no convex hull")

exit()

for buffer in np.arange(0,6,0.5):

    convex_hull_ignore_nearest = False
    params.exhaustive.options.buffer = buffer
    params.exhaustive.output.csv_name = "{}_convex_hull_buffer_ignore_nearest.csv".format(
        str(params.exhaustive.options.buffer).replace('.','_'))
    exhaustive(params=params)