import csv
import os
import logging
import sqlite3

from exhaustive.exhaustive.exhaustive import run as exhaustive
from exhaustive.exhaustive.plotting.plot import scatter_plot
from exhaustive.exhaustive.utils.utils import get_minimum_fofc, u_iso_to_b_fac
from phil import master_phil
from plot_select_regions import plot_protein_region

params =  master_phil.extract()

def list_files(directory, extension):
    return [f for f in os.listdir(directory) if f.endswith('.' + extension)]

def parse_repeat_soak_csv(params):

    input_df = pd.read_csv(params.input.csv)
    for index, row in input_df.iterrows():
        yield row["CrystalName"],row["RefinementPDB_latest"], row["RefinementMTZ_latest"]

def get_xtals_from_db(params,
                      refinement_outcomes="'3 - In Refinement',"
                                          "'4 - CompChem ready', "
                                          "'5 - Deposition ready',"
                                          "'6 - Deposited'"):

    assert os.path.isfile(params.input.database_path), \
        "The database file: \n {} \n does not exist".format(params.input.database_path)

    # Open connection to sqlite database
    conn = sqlite3.connect(params.input.database_path)
    cur = conn.cursor()

    cur.execute("SELECT CrystalName, RefinementPDB_latest, RefinementMTZ_latest "
                "FROM mainTable WHERE RefinementOutcome in ({})" 
                " AND  (RefinementPDB_latest AND RefinementMTZ_latest) IS NOT NULL".format(refinement_outcomes))

    refinement_xtals = cur.fetchall()

    # Close connection to the database
    cur.close()

    for xtal_name, pdb, mtz in refinement_xtals:
        pdb = pdb.encode('ascii')
        mtz = mtz.encode('ascii')
        xtal_name = xtal_name.encode('ascii')
        yield xtal_name, pdb, mtz


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

# out_dir =  "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/FALZA_exhaus_18_09_18_step_0_01_low_U_iso_0/"
# #loop_dir= "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/"
# loop_dir = "/dls/labxchem/data/2016/lb13385-61/processing/analysis/initial_model"
# if not os.path.exists(out_dir):
#     os.mkdir(out_dir)

# xtals=['FALZA-x0079','FALZA-x0085','FALZA-x0172','FALZA-x0177','FALZA-x0271','FALZA-x0309','FALZA-x0402','FALZA-x0438']

# for num in range(start_xtal_num, end_xtal_num + 1):
#     xtal_name = prefix + "{0:0>4}".format(num)
#     xtals.append(xtal_name)
#
# print(xtals)
out_dir =  "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/DCP2B_18_09_20_exhaus"



params = master_phil.extract()

params.settings.processes = 14
params.exhaustive.options.step = 0.01
params.exhaustive.options.convex_hull = False
params.exhaustive.options.per_residue = True
params.exhaustive.options.ligand_grid_points = False
params.exhaustive.options.generate_mtz = False
params.exhaustive.options.lower_u_iso = 0.00
params.input.database_path = "/dls/labxchem/data/2016/lb13385-64/processing/database/soakDBDataFile.sqlite"

if not os.path.exists(out_dir):
    logging.info('Creating output directory {}'.format(out_dir))
    os.mkdir(out_dir)
else:
    logging.info('Output directory {} exists and is being used'.format(params.output.out_dir))

logging.info('Looping over all files that are \'in refinement\' '
            'or better in the supplied datafile: \n {}'.format(params.input.database_path))

for xtal_name, pdb, mtz in get_xtals_from_db(params,
                                             refinement_outcomes="'4 - CompChem ready', "
                                                                 "'5 - Deposition ready',"
                                                                 "'6 - Deposited'" ):

    logging.info(xtal_name)

    assert os.path.exists(pdb), 'PDB File does not exist: {}'.format(pdb)
    assert os.path.exists(mtz), 'MTZ File does not exist: {}'.format(mtz)

    params.input.xtal_name = xtal_name
    params.input.pdb = pdb
    params.input.mtz = mtz
    params.exhaustive.output.csv_name = "exhaustive_search.csv"
    params.output.out_dir = os.path.join(out_dir, xtal_name)

    if not os.path.exists(params.output.out_dir):
        os.mkdir(os.path.join(params.output.out_dir))

    os.chdir(os.path.join(params.output.out_dir))

    try:
        exhaustive(params)
    except UnboundLocalError:
        logging.info("Skipping onto the next crystal")
        continue


    scatter_plot(params.exhaustive.output.csv_name)

    logging.info('Completed: {}'.format(xtal_name))


# xtal_dirs = [os.path.join(loop_dir,xtal_dir) for xtal_dir in os.listdir(loop_dir)
#              if os.path.isdir(os.path.join(loop_dir, xtal_dir))]
#
# csv_paths = []
#
# for xtal_dir in xtal_dirs:
#
#     xtal_name = os.path.basename(xtal_dir)
#     print(xtal_name)
#
#
#     if xtal_name in xtals:
#
#         params.input.xtal_name = xtal_name
#
#         compounds = list_files(xtal_dir,"cif")
#         compound_name = (list_files(xtal_dir,"cif")[0]).split(".")[0]
#
#         params.input.pdb = os.path.join(xtal_dir,"refine.pdb")
#         params.input.mtz = os.path.join(xtal_dir,"refine.mtz")
#         params.output.out_dir = os.path.join(out_dir, compound_name, xtal_name)
#
#         if not os.path.exists(os.path.join(out_dir, compound_name)):
#             os.mkdir(os.path.join(out_dir, compound_name))
#
#         if not os.path.exists(params.output.out_dir):
#             os.mkdir(params.output.out_dir)
#
#         if not os.path.exists(params.input.pdb):
#             print("input pdb doesn't exist: {}".format(params.input.pdb))
#             continue
#         if not os.path.exists(params.input.mtz):
#             print("input mtz doesn't exsit: {}".format(params.input.mtz))
#             continue
#
#         params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")
#
#         csv_paths.append(params.exhaustive.output.csv_name)
#
#         exhaustive(params=params)
#         scatter_plot(params.exhaustive.output.csv_name)



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
#             print([xtal_name, occ, b_fac, fofc])
#
#             minima_writer.writerow([xtal_name, occ, b_fac, fofc])


with open(os.path.join(out_dir,"es_minima.csv"),'wb') as minima_csv:

    minima_writer = csv.writer(minima_csv, delimiter=',')

    for path in csv_paths:
        occ, u_iso, fofc = get_minimum_fofc(path)
        b_fac = u_iso_to_b_fac(u_iso)

        xtal_name = os.path.split(os.path.split(path)[0])[1]

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
  