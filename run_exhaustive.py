import csv
import os

from exhaustive.exhaustive.exhaustive import run as exhaustive
from exhaustive.exhaustive.plotting.plot import scatter_plot
from exhaustive.exhaustive.utils.utils import get_minimum_fofc, u_iso_to_b_fac
from phil import master_phil
from plot_select_regions import plot_protein_region

params =  master_phil.extract()

# example for a single dataset

# params.input.pdb = "/dls/labxchem/data/2018/lb18145-55/processing/analysis/initial_model/NUDT22A-x0927/refine.pdb"
# params.input.mtz = "/dls/labxchem/data/2018/lb18145-55/processing/analysis/initial_model/NUDT22A-x0927/refine.mtz"
# params.input.xtal_name = "NUDT22A-x0927"
# params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
#                         "exhaustive_search_data/test_occ_group_states"
# params.output.log_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
#                         "exhaustive_search_data/test_occ_group_states/logs"
# params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "NUDT22A-x0927_test_occ_group.csv")

# params.input.pdb = "/dls/labxchem/data/2016/lb13385-61/processing/analysis/initial_model/FALZA-x0085/refine.pdb"
# params.input.mtz = "/dls/labxchem/data/2016/lb13385-61/processing/analysis/initial_model/FALZA-x0085/refine.mtz"
# params.input.xtal_name = "FALZA-x0085"
# params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
#                          "exhaustive_search_data/FALZA_exhaus_18_09_14"
# params.output.log_dir = os.path.join(params.output.out_dir, "logs")
# params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")

#Running exhaustive search for covalent ratios/ titration series

# start_xtal_num = 6192
# end_xtal_num = 6251
# #in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios"
# # in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_1812_test"
# in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios_exhaus_sep_18"
# out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios_exhaus_sep_18"
# prefix = "NUDT7A-x"
# qsub = False

#Running exhaustive search for covalent ratios dose experiements

# start_xtal_num = 6192
# end_xtal_num = 6251
# out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios_dose"
# prefix = "NUDT7A-x"
qsub = False

# # copy data to new folder

# if not os.path.exists(out_dir):
#     os.mkdir(out_dir)
#     os.system('cp -a {}/. {}'.format(in_dir,out_dir))

# Single dataset
# params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")
# exhaustive(params=params)
# scatter_plot(params.exhaustive.output.csv_name)
# plot_protein_region(params)

#FALZA exhaustive search

out_dir =  "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/FALZA_exhaus_18_09_18/"
#loop_dir= "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/"
loop_dir = "/dls/labxchem/data/2016/lb13385-61/processing/analysis/initial_model"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

xtals=['FALZA-x0079','FALZA-x0085','FALZA-x0172','FALZA-x0177','FALZA-x0271','FALZA-x0309','FALZA-x0402','FALZA-x0438']

# Multiprocessing using libtbx.easy_mp seems to be failing
params.settings.processes = 14
params.exhaustive.options.step = 0.05
params.exhaustive.options.convex_hull = False
params.exhaustive.options.per_residue = True
params.exhaustive.options.ligand_grid_points = False
params.exhaustive.options.generate_mtz = True


xtals = []
# for num in range(start_xtal_num, end_xtal_num + 1):
#     xtal_name = prefix + "{0:0>4}".format(num)
#     xtals.append(xtal_name)
#
# print(xtals)

xtal_dirs = [os.path.join(loop_dir,xtal_dir) for xtal_dir in os.listdir(loop_dir)
             if os.path.isdir(os.path.join(loop_dir, xtal_dir))]
for xtal_dir in xtal_dirs:

    xtal_name = os.path.basename(xtal_dir)

    if xtal_name in xtals:

        compounds = list_files(xtal_dir,"cif")
        compound_name = (list_files(xtal_dir,"cif")[0]).split(".")[0]

        params.input.pdb = os.path.join(xtal_dir,"refine.pdb")
        params.input.mtz = os.path.join(xtal_dir,"refine.mtz")
        params.output.out_dir = os.path.join(out_dir, compound_name, xtal_name)

        if not os.path.exists(params.input.pdb):
            print("input pdb doesn't exist: {}".format(params.input.pdb))
            continue
        if not os.path.exists(params.input.mtz):
            print("input mtz doesn't exsit: {}".format(params.input.mtz))
            continue

        params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")

        exhaustive(params=params)
        scatter_plot(params.exhaustive.output.csv_name)



# for xtal_name in xtals:
#
#     params.input.xtal_name = xtal_name
#     params.input.pdb = os.path.join(os.path.join(in_dir, xtal_name, "refine.pdb"))
#     params.input.mtz = os.path.join(os.path.join(in_dir, xtal_name, "refine.mtz"))
#
#     if not os.path.exists(params.input.pdb):
#         print("input pdb doesn't exist: {}".format(params.input.pdb))
#         continue
#     if not os.path.exists(params.input.mtz):
#         print("input mtz doesn't exsit: {}".format(params.input.mtz))
#         continue
#
#     params.output.out_dir = os.path.join(out_dir, xtal_name)
#     params.output.log_dir = os.path.join(out_dir, xtal_name, "logs")
#     params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")
#
#     if not qsub:
#         exhaustive(params=params)
#         scatter_plot(params.exhaustive.output.csv_name)
#
#     if qsub:
#
#         # pickle params
#
#         with open(os.path.join(out_dir,
#                                xtal_name,
#                                '{}_param.pickle'.format(xtal_name)),'wb') as handle:
#             pickle.dump(params, handle, protocol=pickle.HIGHEST_PROTOCOL)
#
#         # write python script
#
#         python_script = open(os.path.join(out_dir,
#                                           xtal_name,
#                                           "{}_exhaustive.py".format(xtal_name)),'w')
#         python_script.write("import pickle\n"
#                             "import sys\n"
#                             "import os\n"
#                             "sys.path.append(\"/dls/science/groups/i04-1/elliot-dev/"
#                             "Work/exhaustive_search\")\n"
#                             "from exhaustive.exhaustive import run as exhaustive\n"
#                             "import libtbx.phil\n"
#                             "out_dir=\"{}\"\n".format(out_dir) +
#                             "xtal_name=\"{}\"\n".format(xtal_name) +
#                             "with open(os.path.join(out_dir, " \
#                             "xtal_name,\'{}_param.pickle\'.format(xtal_name)),'rb') as handle:\n"
#                             "\tparams = pickle.load(handle)\n"
#                             "exhaustive(params)")
#         python_script.close()
#
#         # write bash script
#         bash_script=open(os.path.join(out_dir, xtal_name,
#                                      "{}_exhaustive.sh".format(xtal_name)),'w')
#         bash_script.write("source /dls/science/groups/i04-1/software/" \
#                           "pandda-update/ccp4/ccp4-7.0/setup-scripts/ccp4.setup-sh\n"
#                           # "export PYTHONPATH=\"${PYTHONPATH}:" \
#                           # "/dls/science/groups/i04-1/elliot-dev/"
#                           # "Work/exhaustive_search/exhaustive\"\n"
#                           "ccp4-python " + os.path.join(out_dir, xtal_name,
#                                           "{}_exhaustive.py".format(xtal_name)))
#         bash_script.close()
#         # submit job
#         os.system("qsub {}".format(os.path.join(out_dir, xtal_name,xtal_name+"_exhaustive.sh")))


# Get exhaustive search minima fofc
with open(os.path.join(out_dir,"es_minima.csv"),'wb') as minima_csv:

    minima_writer = csv.writer(minima_csv, delimiter=',')

    for xtal_name in xtals:

        params.output.out_dir = os.path.join(out_dir, xtal_name)
        params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")
        if os.path.exists(params.exhaustive.output.csv_name):
            os.chdir(os.path.join(out_dir, xtal_name))
            scatter_plot(params.exhaustive.output.csv_name)
        else:
            continue

        if os.path.exists(params.exhaustive.output.csv_name):
            occ, u_iso, fofc = get_minimum_fofc(params.exhaustive.output.csv_name)
            b_fac=u_iso_to_b_fac(u_iso)

            print([xtal_name, occ, b_fac, fofc])

            minima_writer.writerow([xtal_name, occ, b_fac, fofc])


#refine minima

# with open(os.path.join(out_dir,"refined_occs.csv"),'wb') as minima_csv:
#
#     minima_writer = csv.writer(minima_csv, delimiter=',')
#
#     for xtal_name in xtals:
#
#         if os.path.exists(os.path.join(out_dir,xtal_name,"refine.pdb")):
#
#             occ = get_lig_occ(os.path.join(out_dir,xtal_name,"refine.pdb"))
#
#             minima_writer.writerow([xtal_name,occ])
#         else:
#             continue
  