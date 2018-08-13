from exhaustive import run as exhaustive
from phil import master_phil
from plotting.plot import scatter_plot
import numpy as np
import os

params =  master_phil.extract()
params.input.xtal_name = "NUDT22A-x"

params.input.pdb = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/" \
"repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/" \
"FMOPL000622a_DSI_poised/NUDT22A-x0955/refine.pdb"
params.input.mtz = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/" \
"repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/" \
"FMOPL000622a_DSI_poised/NUDT22A-x0955/refine.mtz"

params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/convex_buffer_tests/NUDT22A-x0955"
params.settings.processes = 1

# params.exhaustive.output.csv_name = "no_convex_hull.csv"
# params.exhaustive.options.convex_hull = False
# exhaustive(params = params)

scatter_plot(os.path.join(params.output.out_dir,params.exhaustive.output.csv_name), three_dim_plot=True)

params.exhaustive.options.convex_hull = True
params.exhaustive.output.csv_name = "convex_hull.csv"
exhaustive(params = params)
scatter_plot(os.path.join(params.output.out_dir,params.exhaustive.output.csv_name), three_dim_plot=True)
exit()

exit()


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