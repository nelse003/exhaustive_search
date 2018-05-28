import os
import datetime
import sqlite3

from exhaustive.exhaustive import run as exhaustive
from exhaustive.phil import master_phil, check_input_files


params =  master_phil.extract()

# Setting common paramters
common_out_dir = os.path.join("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/",
                        datetime.datetime.now().strftime('%Y-%m-%d')) #_%H-%M-%S

if not os.path.exists(common_out_dir):
    os.mkdir(common_out_dir)

params.exhaustive.options.generate_mtz = False
params.exhaustive.options.step = 0.01
params.settings.processes = 36

# Running exhaustive search

# NUDT7 copied atoms (OX210)

# NUDT7_copied_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/NUDT7_Copied_atoms"
#
# for NUDT7_cov_dataset in filter(lambda x: x.startswith("NUDT7A-x") and
#                        os.path.isdir(os.path.join(NUDT7_copied_dir, x)), os.listdir(NUDT7_copied_dir)):
#
#     params.input.xtal_name = NUDT7_cov_dataset
#     params.output.out_dir = os.path.join(common_out_dir, "NUDT7_copied", params.input.xtal_name)
#
#     if not os.path.exists(os.path.join(common_out_dir, "NUDT7_copied")):
#         os.mkdir(os.path.join(common_out_dir, "NUDT7_copied"))
#
#     if not os.path.exists(os.path.join(common_out_dir, "NUDT7_copied",params.input.xtal_name)):
#         os.mkdir(os.path.join(common_out_dir, "NUDT7_copied",params.input.xtal_name))
#
#     params.input.in_path = os.path.join(NUDT7_copied_dir, NUDT7_cov_dataset)
#     params.input.mtz = os.path.join(params.input.in_path, "refine.mtz")
#     params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")
#     params.output.log_dir = os.path.join(params.output.out_dir, "logs")
#
#     params.exhaustive.output.csv_name = params.input.xtal_name + "exhaustive_search_occ_u_iso"
#
#     if not os.path.exists(params.output.log_dir):
#         os.mkdir(params.output.log_dir)
#
#     print(params.input.xtal_name)
#     check_input_files(params)
#     exhaustive(params)

# NUDT7 NUOOA0000181a

#
# # NUDT 7 Covalent hits?
#
# # NUDT22A FMOPL00622a
# from DSI poised library (x0938-x0976)
# from DSPL library x0182, x909-x0937

# for in:
#     params.output.out_dir =
#     params.input.xtal_name =
#     params.input.in_path =
#     params.input.mtz =
#     params.input.pdb =
#     params.output.log_dir = os.path.join(params.output.out_dir, "logs")
#
#     exhaustive(params)

# NUDT22A 133725a x0421, x1040 - x1059
# NUDT22A 13663a x0391, x1009 to x1039
# NUDT22A 13369a x0243, x977 to x1008

NUDT22_initial_models = "/dls/labxchem/data/2018/lb18145-55/processing/analysis/initial_model"
params.repeat.input.database_path = "/dls/labxchem/data/2018/lb18145-55/processing/database/soakDBDataFile.sqlite"

rejects = []

start_xtal_num = 977
end_xtal_num = 1058
prefix = "NUDT22A-x"
xtals = ['NUDT22A-x0243', 'NUDT22A-x0421','NUDT22A-x0391']
xtals =[]
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = prefix + "{0:0>4}".format(num)
    xtals.append(xtal_name)

for xtal_name in xtals:

    params.input.xtal_name = xtal_name
    params.input.in_path = os.path.join(NUDT22_initial_models,xtal_name)
    params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")
    print(params.input.pdb)
    params.output.out_dir = os.path.join(common_out_dir, "NUDT22_from_NUDT22A_Elliot", xtal_name)

    if not os.path.exists(os.path.join(common_out_dir, "NUDT22_from_NUDT22A_Elliot")):
        os.mkdir(os.path.join(common_out_dir, "NUDT22_from_NUDT22A_Elliot"))

    if not os.path.exists(os.path.join(common_out_dir, "NUDT22_from_NUDT22A_Elliot",params.input.xtal_name)):
        os.mkdir(os.path.join(common_out_dir, "NUDT22_from_NUDT22A_Elliot",params.input.xtal_name))

    params.output.log_dir = os.path.join(params.output.out_dir, "logs")
    if not os.path.exists(params.output.log_dir):
        os.mkdir(params.output.log_dir)

    params.exhaustive.output.csv_name = params.input.xtal_name + "exhaustive_search_occ_u_iso"

    if os.path.exists(os.path.join(params.output.out_dir,params.exhaustive.output.csv_name)):
        continue


    # Get mtz from database
    conn = sqlite3.connect(params.repeat.input.database_path)
    cur = conn.cursor()
    cur.execute("SELECT RefinementMTZ_latest "
                "FROM mainTable WHERE CrystalName=?",(xtal_name,) )

    refinement_xtals = cur.fetchall()

    # Close connection to the database
    cur.close()

    if refinement_xtals[0][0] is not None:
        params.input.mtz = refinement_xtals[0][0].encode('ascii')

        check_input_files(params)
        try:
            exhaustive(params)
        except UnboundLocalError:
            rejects.append(xtal_name)
            continue

    else:
        print("Refinement Mtz does not exist? for xtal: {}".format(xtal_name))
        rejects.append(xtal_name)
        continue

print(rejects)

# DCP2B FMOPL000435a
#
# for in:
#
#     params.output.out_dir =
#     params.input.xtal_name =
#     params.input.in_path =
#     params.input.mtz =
#     params.input.pdb =
#     params.output.log_dir = os.path.join(params.output.out_dir, "logs")
#
#     exhaustive(params)

