from exhaustive import run as exhaustive
from phil import master_phil
from plotting.plot import scatter_plot
import os
import pickle
import csv
from utils.utils import get_minimum_fofc, u_iso_to_b_fac
params =  master_phil.extract()

# example for a single dataset

# params.input.pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/" \
#                    "initial_model/NUDT7A-x0299/refine.pdb"
# params.input.mtz = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/" \
#                    "initial_model/NUDT7A-x0299/refine.mtz"
# params.input.xtal_name = "NUDT7A-x0299"
# params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
#                         "exhaustive_search_data/convex_hull"
# params.output.log_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
#                         "exhaustive_search_data/convex_hull/logs"
# params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "result_convex_buffer.csv")

# Multiprocessing using libtbx.easy_mp seems to be failing
params.settings.processes = 1
params.exhaustive.options.step = 0.01
params.exhaustive.options.convex_hull = False


#Running exhaustive search for covalent ratios/ titration series

start_xtal_num = 1905
end_xtal_num = 2005
out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratio_no_convex"
prefix = "NUDT7A-x"
qsub = True

xtals = []
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = prefix + "{0:0>4}".format(num)
    xtals.append(xtal_name)

for xtal_name in xtals:

    params.input.xtal_name = xtal_name
    params.input.pdb = os.path.join(os.path.join(out_dir, xtal_name, "refine.pdb"))
    params.input.mtz = os.path.join(os.path.join(out_dir, xtal_name, "refine.mtz"))

    if not os.path.exists(params.input.pdb):
        continue
    if not os.path.exists(params.input.mtz):
        continue

    params.output.out_dir = os.path.join(out_dir, xtal_name)
    params.output.log_dir = os.path.join(out_dir, xtal_name, "logs")
    params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")

    if not qsub:
        exhaustive(params=params)
        scatter_plot(params.exhaustive.output.csv_name)

    if qsub:

        # pickle params

        with open(os.path.join(out_dir,
                               xtal_name,
                               '{}_param.pickle'.format(xtal_name)),'wb') as handle:
            pickle.dump(params, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # write python script

        python_script = open(os.path.join(out_dir,
                                          xtal_name,
                                          "{}_exhaustive.py".format(xtal_name)),'w')
        python_script.write("import pickle\n"
                            "import sys\n"
                            "import os\n"
                            "sys.path.append(\"/dls/science/groups/i04-1/elliot-dev/"
                            "Work/exhaustive_search\")\n"
                            "from exhaustive.exhaustive import run as exhaustive\n"
                            "import libtbx.phil\n"
                            "out_dir=\"{}\"\n".format(out_dir) +
                            "xtal_name=\"{}\"\n".format(xtal_name) +
                            "with open(os.path.join(out_dir, " \
                            "xtal_name,\'{}_param.pickle\'.format(xtal_name)),'rb') as handle:\n"
                            "\tparams = pickle.load(handle)\n"
                            "exhaustive(params)")
        python_script.close()

        # write bash script
        bash_script=open(os.path.join(out_dir, xtal_name,
                                     "{}_exhaustive.sh".format(xtal_name)),'w')
        bash_script.write("source /dls/science/groups/i04-1/software/" \
                          "pandda-update/ccp4/ccp4-7.0/setup-scripts/ccp4.setup-sh\n"
                          # "export PYTHONPATH=\"${PYTHONPATH}:" \
                          # "/dls/science/groups/i04-1/elliot-dev/"
                          # "Work/exhaustive_search/exhaustive\"\n"
                          "ccp4-python " + os.path.join(out_dir, xtal_name,
                                          "{}_exhaustive.py".format(xtal_name)))
        bash_script.close()
        # submit job
        os.system("qsub {}".format(os.path.join(out_dir, xtal_name,xtal_name+"_exhaustive.sh")))

# Get exhaustive search minima fofc
# with open(os.path.join(out_dir,"es_minima.csv"),'wb') as minima_csv:
#
#     minima_writer = csv.writer(minima_csv, delimiter=',')
#
#     for xtal_name in xtals:
#
#         params.output.out_dir = os.path.join(out_dir, xtal_name)
#         params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")
#         if os.path.exists(params.exhaustive.output.csv_name):
#             os.chdir(os.path.join(out_dir, xtal_name))
#             scatter_plot(params.exhaustive.output.csv_name)
#         else:
#             continue
#
#         if os.path.exists(params.exhaustive.output.csv_name):
#             occ, u_iso, fofc = get_minimum_fofc(params.exhaustive.output.csv_name)
#             b_fac=u_iso_to_b_fac(u_iso)
#
#             minima_writer.writerow([xtal_name, occ, b_fac, fofc])



  