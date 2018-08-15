from exhaustive import run as exhaustive
from phil import master_phil
from plotting.plot import scatter_plot
import numpy as np
import os

params =  master_phil.extract()
params.input.xtal_name = "NUDT22A-x0955"

# params.input.pdb = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/" \
# "NUDT7_covalent/NUDT7A-x1812/refine_0007/output.pdb"
# params.input.mtz = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/" \
# "NUDT7_covalent/NUDT7A-x1812/refine_0007/output.mtz"

params.input.pdb = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_search_data/repeat_soaks/2018-05-28/" \
                   "NUDT22_from_occ_group_with_refinement/" \
                   "FMOPL000622a_DSI_poised/NUDT22A-x0955/refine.pdb"
params.input.mtz = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_search_data/repeat_soaks/2018-05-28/" \
                   "NUDT22_from_occ_group_with_refinement/" \
                   "FMOPL000622a_DSI_poised/NUDT22A-x0955/refine.mtz"

params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                        "exhaustive_search_data/convex_buffer_tests/" \
                        "NUDT22A-x0955-after-num-altlocs/mtz-test"

params.settings.processes = 1

# params.exhaustive.output.csv_name = "no_convex_hull.csv"
# params.exhaustive.options.convex_hull = False
# exhaustive(params = params)
# scatter_plot(os.path.join(params.output.out_dir,params.exhaustive.output.csv_name), three_dim_plot=True)

params.exhaustive.options.generate_mtz= True
params.exhaustive.options.convex_hull = True
params.exhaustive.output.csv_name = "convex_hull.csv"
exhaustive(params = params)
scatter_plot(os.path.join(params.output.out_dir,params.exhaustive.output.csv_name), three_dim_plot=True)

# params.exhaustive.options.convex_hull = False
# params.exhaustive.options.ligand_atom_points = True
# params.exhaustive.options.atom_points_sel_string = \
#     "(chain B and altid C and resid 1) or (chain B and altid D resid 1)"
# params.exhaustive.output.csv_name = "lig_atom_points_convex_hull.csv"
# exhaustive(params = params)
# scatter_plot(os.path.join(params.output.out_dir,params.exhaustive.output.csv_name), three_dim_plot=True)
#
# params.exhaustive.options.convex_hull = False
# params.exhaustive.options.ligand_grid_points = True
# params.exhaustive.options.atom_points_sel_string = \
#     "(chain B and altid C and resid 1) or (chain B and altid D resid 1)"
# params.exhaustive.output.csv_name = "lig_grid_points_convex_hull.csv"
# exhaustive(params = params)
# scatter_plot(os.path.join(params.output.out_dir,params.exhaustive.output.csv_name), three_dim_plot=True)
#
# for buffer in np.arange(0,4,0.5):

#     params.exhaustive.options.convex_hull = True
#     params.exhaustive.options.buffer = buffer
#     params.exhaustive.output.csv_name = "{}_convex_hull_buffer.csv".format(
#         str(params.exhaustive.options.buffer).replace('.','_'))
#     exhaustive(params=params)
#     scatter_plot(os.path.join(params.output.out_dir,
#                           params.exhaustive.output.csv_name),
#                               three_dim_plot=True,
#                  title_text="{}_convex_hull_buffer".format(str(
#                      params.exhaustive.options.buffer).replace('.', '_')))
#
# scatter_plot(os.path.join(params.output.out_dir,
#                           params.exhaustive.output.csv_name),
#              three_dim_plot=True,
#              title_text="convex_hull_no_buffer")
#
# for buffer in np.arange(0,4,0.5):
#
#     convex_hull_ignore_nearest = True
#     params.exhaustive.options.buffer = buffer
#     params.exhaustive.output.csv_name = "{}_convex_hull_buffer_ignore_nearest.csv".format(
#         str(params.exhaustive.options.buffer).replace('.','_'))
#     exhaustive(params=params)
#     scatter_plot(os.path.join(params.output.out_dir,
#                               params.exhaustive.output.csv_name),
#                  three_dim_plot=True,
#                  title_text="convex hull ignore nearest atoms buffer {}".format(params.exhaustive.options.buffer))