from __future__ import division, print_function
import os
import datetime
import sqlite3
import shutil
import itertools
import libtbx.phil

from exhaustive.exhaustive import run as exhaustive
from exhaustive.phil import master_phil, check_input_files
from exhaustive.utils.utils import get_minimum_fofc, get_occ_b
from exhaustive.utils.utils import u_iso_to_b_fac
from exhaustive.plotting.plot import occupancy_histogram_with_exhaustive_search
from exhaustive.plotting.plot import occupancy_b_factor_scatter_plot
from exhaustive.process.minima import write_minima_pdb

from giant.jiffies.split_conformations import run as split_conformations
from giant.jiffies.split_conformations import master_phil as split_phil
from giant.jiffies.merge_conformations import run as merge_conformations
from giant.jiffies.merge_conformations import master_phil as merge_phil
from giant.jiffies.score_model import run as score_model
from giant.jiffies.score_model import master_phil as score_phil
params =  master_phil.extract()

#TODO Apply DRY and functionalise this script

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

#######################################################
# Run exhaustive search @ 0.01
#######################################################

NUDT7_copied_dir = ("/dls/science/groups/i04-1/elliot-dev/Work/"
                   "exhaustive_search_data/NUDT7_Copied_atoms")

# for NUDT7_dataset in filter(lambda x: x.startswith("NUDT7A-x") and
#                        os.path.isdir(os.path.join(NUDT7_copied_dir, x)), os.listdir(NUDT7_copied_dir)):
#
#    params.input.xtal_name = NUDT7_dataset
#    params.output.out_dir = os.path.join("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-27"
#                                         , "NUDT7_copied", params.input.xtal_name)
#
#     if not os.path.exists(os.path.join(common_out_dir, "NUDT7_copied")):
#         os.mkdir(os.path.join(common_out_dir, "NUDT7_copied"))
#
#     if not os.path.exists(os.path.join(common_out_dir, "NUDT7_copied",params.input.xtal_name)):
#         os.mkdir(os.path.join(common_out_dir, "NUDT7_copied",params.input.xtal_name))
#
#     params.input.in_path = os.path.join(NUDT7_copied_dir, NUDT7_dataset)
#     params.input.mtz = os.path.join(params.input.in_path, "refine.mtz")
#     params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")
#     params.output.log_dir = os.path.join(params.output.out_dir, "logs")
#     params.exhaustive.output.csv_name = params.input.xtal_name + "_exhaustive_search_occ_u_iso.csv"
#
#     if not os.path.exists(params.output.log_dir):
#         os.mkdir(params.output.log_dir)
#
#     print(params.input.xtal_name)
#     check_input_files(params)
#     exhaustive(params)

################################################################
# Plot the copied atoms refinement and exhaustive search occupancy
################################################################

# refine_occs, mean_ligand_b_factor, std_ligand_b_fac = get_occ_b(
#     refinement_dir=NUDT7_copied_dir,
#     lig_chain="E",
#     pdb_name="refine.split.bound-state_with_new_atoms.pdb")
#
# # Adding in extra dataset
# refine_occs_extra, \
# mean_ligand_b_factor_extra, \
# std_ligand_b_fac_extra = get_occ_b(refinement_dir=NUDT7_copied_dir,
#                                     lig_chain="E",
#                                     pdb_name="refine.split.bound-state.pdb")
# refine_occs = refine_occs + refine_occs_extra
# mean_ligand_b_factor = mean_ligand_b_factor + mean_ligand_b_factor_extra
# std_ligand_b_fac = std_ligand_b_fac + std_ligand_b_fac_extra
#
# es_occs =[]
# es_b_fac = []
# for NUDT7_cov_dataset in filter(lambda x: x.startswith("NUDT7A-x") and
#                        os.path.isdir(
#                            os.path.join("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-27/NUDT7_copied", x)),
#                             os.listdir("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-27/NUDT7_copied")):
#
#     csv_name = os.path.join("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-27/NUDT7_copied",
#                             NUDT7_cov_dataset, NUDT7_cov_dataset + "_exhaustive_search_occ_u_iso.csv")
#     print(csv_name)
#     occ, u_iso , _ = get_minimum_fofc(csv_name)
#     es_b_fac.append(u_iso_to_b_fac(u_iso))
#     es_occs.append(occ)
#
#
# params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-27/NUDT7_copied"
# #occupancy_histogram_with_exhaustive_search(es_occs, refine_occs,protein_name="NUDT7", compound="OX210", params=params)
#
#
# occupancy_b_factor_scatter_plot(es_occs = es_occs,
#                                 refined_occs = refine_occs,
#                                 es_b_fac = es_b_fac ,
#                                 refine_mean_b_fac = mean_ligand_b_factor,
#                                 refine_std_b_fac = std_ligand_b_fac,
#                                 protein_name = "NUDT7A",
#                                 compound = "OX210",
#                                 params = params)

# NUDT7 Covalent

# Going to run on the files that are in
# /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/NUDT7_covalent
# Uncertain whether these come from copied atoms, the multi-state-model-pdb's
# are different, so we assume these come from copied files which
# are then refined. These refinements are in refine_0006


# NUDT7_covalent_path = "/dls/science/groups/i04-1/elliot-dev/" \
#                      "Work/exhaustive_search_data/NUDT7_covalent"
#
# if not os.path.exists(os.path.join(common_out_dir, "NUDT7_cov")):
#     os.mkdir(os.path.join(common_out_dir, "NUDT7_cov"))
#
# for NUDT7_cov_dataset in filter(lambda x: x.startswith("NUDT7A-x") and
#                        os.path.isdir(os.path.join(NUDT7_covalent_path, x)),
#                                 os.listdir(NUDT7_covalent_path)):
#
#     # Get location of most recent output.pdb
#     output_pdb_folders = []
#     output_pdb_path = None
#     for folder in filter(lambda x: os.path.isdir(
#         os.path.join(NUDT7_covalent_path,x)),
#         os.listdir(os.path.join(NUDT7_covalent_path,NUDT7_cov_dataset))):
#
#         for file in filter(lambda x: x == "output.pdb" and os.path.isfile(
#                 os.path.join(NUDT7_covalent_path, NUDT7_cov_dataset,
#                              folder, x)),
#                            os.listdir(os.path.join(NUDT7_covalent_path,
#                                                    NUDT7_cov_dataset, folder))):
#
#             output_pdb_folders.append(int(folder[-4:]))
#
#     if output_pdb_folders:
#
#         folder = "refine_" + str(max(output_pdb_folders)).zfill(4)
#
#         output_pdb_path = os.path.join(NUDT7_covalent_path,
#                                        NUDT7_cov_dataset,
#                                        folder,
#                                        "output.pdb")
#
#         print("PATH: " + output_pdb_path)
#
#     if not os.path.exists(os.path.join(NUDT7_covalent_path, NUDT7_cov_dataset,
#                                        "refine.split.ground-state.pdb")):
#         merge_params = merge_phil.extract()
#
#         merge_params.input.minor = os.path.join(
#             "/dls/labxchem/data/2017/lb18145-49/processing/"
#             "analysis/initial_model/", NUDT7_cov_dataset,
#             NUDT7_cov_dataset +"-pandda-model.pdb")
#
#         merge_params.input.major = os.path.join(
#             "/dls/labxchem/data/2017/lb18145-49/processing/"
#             "analysis/initial_model/", NUDT7_cov_dataset,
#             NUDT7_cov_dataset +"-pandda-input.pdb")
#
#         merge_params.output.pdb = os.path.join(NUDT7_covalent_path,
#                                                NUDT7_cov_dataset,
#                                                "merged.pdb")
#
#         if not os.path.exists(merge_params.input.major):
#             print("Skipping {} as {} doesn't exist".format(NUDT7_cov_dataset,
#                                                     merge_params.input.major))
#             continue
#
#         if not os.path.exists(merge_params.input.minor):
#             print("Skipping {} as {} doesn't exist".format(NUDT7_cov_dataset,
#                                                     merge_params.input.minor))
#             continue
#
#         print(merge_params.input.minor)
#         print(merge_params.input.major)
#         print(merge_phil.format(python_object=merge_params).as_str())
#         try:
#             merge_conformations(merge_params)
#         except IOError:
#             print(" confs: Issue in parsing:{}".format(NUDT7_cov_dataset))
#             continue
#     exit()
#
#     if output_pdb_path is None:
#         print(os.path.join(NUDT7_covalent_path,NUDT7_cov_dataset))
#         os.chdir(os.path.join(NUDT7_covalent_path,NUDT7_cov_dataset))
#
#         input_pdb = "merged.pdb"
#         input_mtz = ".free.mtz"
#         input_cif = "LIG-CYS.cif"
#
#         try:
#             os.system("giant.quick_refine input.pdb={} "
#                       "input.mtz={} input.cif={}".format(input_pdb,
#                                                          input_mtz,
#                                                          input_cif))
#         except IOError:
#             print "Skipping crystal"
#             continue

# Run exhaustive search

# for NUDT7_cov_dataset in filter(lambda x: x.startswith("NUDT7A-x") and
#                        os.path.isdir(os.path.join(NUDT7_covalent_path, x)),
#                                 os.listdir(NUDT7_covalent_path)):
#
#     params.input.xtal_name = NUDT7_cov_dataset
#     params.output.out_dir = os.path.join(common_out_dir,
#                                         "NUDT7_cov",
#                                         params.input.xtal_name)
#
#     if not os.path.exists(params.output.out_dir):
#         os.mkdir(params.output.out_dir)
#
#     params.input.in_path = os.path.join(NUDT7_covalent_path, NUDT7_cov_dataset)
#     params.input.mtz = os.path.join(params.input.in_path, "refine.mtz")
#     params.input.pdb = output_pdb_path
#     params.output.log_dir = os.path.join(params.output.out_dir, "logs")
#     params.exhaustive.output.csv_name = params.input.xtal_name \
#                                         + "_exhaustive_search_occ_u_iso.csv"
#
#     if not os.path.exists(params.output.log_dir):
#         os.mkdir(params.output.log_dir)
#
#     try:
#         check_input_files(params)
#     except AssertionError:
#         print("Skipping {} as input files aren't ok".format(NUDT7_cov_dataset))
#         continue
#     exhaustive(params)
#
# exit()
################################################################
# NUDT7 Covalent
#
# Plotting occupancy histogram
# Plot the copied atoms refinement and exhaustive search occupancy
################################################################

# refine_occs = get_occs(refinement_dir=NUDT7_covalent_path,lig_chain="E",
#                        pdb_name="output.pdb")
# refine_occs, mean_ligand_b_factor, std_ligand_b_fac = get_occ_b(
#     refinement_dir=NUDT7_covalent_path,
#     lig_chain="E",
#     pdb_name="output.pdb")
#
# NUDT7_covalent_out_path = "/dls/science/groups/i04-1/elliot-dev/Work/" \
#                           "exhaustive_search_data/repeat_soaks/" \
#                           "2018-06-05/NUDT7_cov"
#
# # refine_occs = refine_occs + get_occs(refinement_dir=NUDT7_copied_dir,
# #                           lig_chain="E",
# #                        pdb_name="refine.split.bound-state.pdb")
# #
# es_occs = []
# es_b_fac = []
# for NUDT7_cov_dataset in filter(lambda x: x.startswith("NUDT7A-x") and
#                        os.path.isdir(
#                            os.path.join(NUDT7_covalent_out_path, x)),
#                             os.listdir(NUDT7_covalent_out_path)):
#
#     if not os.path.exists(os.path.join(NUDT7_covalent_path, NUDT7_cov_dataset,
#                                        "refine.split.ground-state.pdb")):
#         split_params = split_phil.extract()
#         split_params.input.pdb = [os.path.join(NUDT7_covalent_path,
#                                                NUDT7_cov_dataset,
#                                                "refine.pdb")]
#         split_params.output.log = "cat.log"
#         print(split_params.input.pdb)
#         print(split_phil.format(python_object=split_params).as_str())
#         try:
#             split_conformations(split_params)
#         except IOError:
#             print("Split confs: Issue in parsing:{}".format(NUDT7_cov_dataset))
#             continue
#
#
#     csv_name = os.path.join("/dls/science/groups/i04-1/elliot-dev/Work/"
#                             "exhaustive_search_data/repeat_soaks/"
#                             "2018-05-27/NUDT7_copied",
#                             NUDT7_cov_dataset,
#                             NUDT7_cov_dataset
#                             + "_exhaustive_search_occ_u_iso.csv")
#     print(csv_name)
#     occ, u_iso, _ = get_minimum_fofc(csv_name)
#     es_b_fac.append(u_iso_to_b_fac(u_iso))
#     es_occs.append(occ)
#
# print(refine_occs,es_occs)
#
# params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
#                         "exhaustive_search_data/repeat_soaks/" \
#                         "2018-05-27/NUDT7_copied"
#
# occupancy_histogram_with_exhaustive_search(es_occs, refine_occs,
#                                            protein_name="NUDT7",
#                                            compound="NUOOA000181a",
#                                            params=params)
#
# occupancy_b_factor_scatter_plot(es_occs = es_occs,
#                                 refined_occs = refine_occs,
#                                 es_b_fac = es_b_fac ,
#                                 refine_mean_b_fac = mean_ligand_b_factor,
#                                 refine_std_b_fac = std_ligand_b_fac,
#                                 protein_name = "NUDT7A",
#                                 compound = "NUOOA000181a",
#                                 params = params)

# NUDT22_initial_models = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/occupancy_group_with_refinement" #"/dls/labxchem/data/2018/lb18145-55/processing/analysis/initial_model"
# params.repeat.input.database_path = "/dls/labxchem/data/2018/lb18145-55/processing/database/soakDBDataFile.sqlite"
#
# rejects = []
#
# start_xtal_num = 909
# end_xtal_num = 1058
# prefix = "NUDT22A-x"
# xtals = ['NUDT22A-x0182','NUDT22A-x0243', 'NUDT22A-x0421','NUDT22A-x0391']
# for num in range(start_xtal_num, end_xtal_num + 1):
#     xtal_name = prefix + "{0:0>4}".format(num)
#     xtals.append(xtal_name)
#
# for xtal_name in xtals:
#
#     params.input.xtal_name = xtal_name
#     params.input.in_path = os.path.join(NUDT22_initial_models,xtal_name)
#     params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")
#     print(params.input.pdb)
#     params.output.out_dir = os.path.join(common_out_dir, "NUDT22_from_occ_group_with_refinement", xtal_name)
#
#     if not os.path.exists(os.path.join(common_out_dir, "NUDT22_from_occ_group_with_refinement")):
#         os.mkdir(os.path.join(common_out_dir, "NUDT22_from_occ_group_with_refinement"))
#
#     if not os.path.exists(os.path.join(common_out_dir,
#                           "NUDT22_from_occ_group_with_refinement",
#                           params.input.xtal_name)):

#         os.mkdir(os.path.join(common_out_dir,
#                               "NUDT22_from_occ_group_with_refinement",
#                               params.input.xtal_name))
#
#     params.output.log_dir = os.path.join(params.output.out_dir, "logs")
#     if not os.path.exists(params.output.log_dir):
#         os.mkdir(params.output.log_dir)
#
#     params.exhaustive.output.csv_name = params.input.xtal_name
#                                         + "exhaustive_search_occ_u_iso"
#
#     if os.path.exists(os.path.join(params.output.out_dir,
#                       params.exhaustive.output.csv_name)):
#         continue
#
#     params.input.mtz = os.path.join(params.input.in_path,"refine.mtz")
#
#     try:
#         check_input_files(params)
#     except:
#         rejects.append(xtal_name)
#         continue
#     # # Get mtz from database
#     # conn = sqlite3.connect(params.repeat.input.database_path)
#     # cur = conn.cursor()
#     # cur.execute("SELECT RefinementMTZ_latest "
#     #             "FROM mainTable WHERE CrystalName=?",(xtal_name,) )
#     #
#     # refinement_xtals = cur.fetchall()
#     #
#     # # Close connection to the database
#     # cur.close()
#     #
#     # if refinement_xtals[0][0] is not None:
#     #     params.input.mtz = refinement_xtals[0][0].encode('ascii')
#     #
#     #     check_input_files(params)
#
#     try:
#         exhaustive(params)
#     except UnboundLocalError:
#         rejects.append(xtal_name)
#         continue
#
# #print(rejects)

######################################################
# NUDT22 Seperate into compound sets
######################################################
prefix = "NUDT22A-x"
xtals = {}
# # NUDT22A FMOPL00622a
# from DSI poised library (x0938-x0976)
compound_code = 'FMOPL000622a_DSPL'
xtals['NUDT22A-x0182'] =  compound_code
start_xtal_num = 909
end_xtal_num = 937
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = prefix + "{0:0>4}".format(num)
    xtals[xtal_name] = compound_code

# # NUDT22A FMOPL00622a
# from DSI poised library (x0938-x0976)
compound_code = 'FMOPL000622a_DSI_poised'
start_xtal_num = 938
end_xtal_num = 976
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = prefix + "{0:0>4}".format(num)
    xtals[xtal_name] = compound_code

# NUDT22A 133725a x0421, x1040 - x1059

start_xtal_num = 1040
end_xtal_num = 1059
compound_code = '133725a'
xtals['NUDT22A-x0421'] = compound_code
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = prefix + "{0:0>4}".format(num)
    xtals[xtal_name] = compound_code

# NUDT22A 13369a x0243, x977 to x1008

start_xtal_num = 977
end_xtal_num = 1008
compound_code = '13369a'
xtals['NUDT22A-x0243'] = compound_code
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = prefix + "{0:0>4}".format(num)
    xtals[xtal_name] = compound_code

# NUDT22A 13663a x0391, x1009 to x1039

compound_code = '13363a'
xtals['NUDT22A-x0391']= compound_code
start_xtal_num = 1009
end_xtal_num = 1039
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = prefix + "{0:0>4}".format(num)
    xtals[xtal_name] = compound_code

in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/occupancy_group_with_refinement"

compound_codes=set()
for compounds in xtals.values():
    compound_codes.add(compounds)

for compound in compound_codes:
    if not os.path.exists(os.path.join(in_dir,compound)):
        os.mkdir(os.path.join(in_dir,compound))

for xtal, compound_code in xtals.iteritems():
    if os.path.exists(os.path.join(in_dir,xtal)):
        shutil.move(os.path.join(in_dir,xtal),os.path.join(in_dir,compound_code))

######################################################
# NUDT22
######################################################

out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/" \
          "repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/"

for compound in compound_codes:
    for dataset in filter(lambda x: x.startswith("NUDT22A-x") and
                           os.path.isdir(os.path.join(in_dir, compound, x)),
                                os.listdir(os.path.join(in_dir, compound))):

        # Generate split conformations
        if not os.path.exists(os.path.join(in_dir, compound, dataset,
                                           "refine.split.ground-state.pdb")):
            split_params = split_phil.extract()
            split_params.input.pdb = [os.path.join(in_dir, compound, dataset,
                                                   "refine.pdb")]
            split_params.output.log = "cat.log"
            print(split_params.input.pdb)
            print(split_phil.format(python_object=split_params).as_str())
            try:
                split_conformations(split_params)
            except IOError:
                print("Split confs: Issue in parsing:{}".format(dataset))
                continue


    refine_occs, mean_ligand_b_factor, std_ligand_b_fac  = get_occ_b(
        refinement_dir=os.path.join(in_dir, compound),
        lig_chain="B",
        pdb_name="refine.split.bound-state.pdb")

    es_occs = []
    es_b_fac = []

    for dataset in filter(lambda x: x.startswith("NUDT22A-x") and os.path.isdir(os.path.join(out_dir, compound, x)),
                          os.listdir(os.path.join(out_dir, compound))):

        csv_name = os.path.join(out_dir, compound, dataset,
                                dataset + "_exhaustive_search_occ_u_iso.csv")

        refine_pdb = os.path.join(in_dir, compound, dataset, "refine.pdb")
        es_pdb = os.path.join(out_dir, compound, dataset, "es_minima.pdb")

        print(refine_pdb)

        # Minima pdb
        if not os.path.exists(es_pdb):
            try:
                write_minima_pdb(input_pdb=refine_pdb,
                                 output_pdb=es_pdb,
                                 csv_name=csv_name,
                                 params=params)
            except IOError:
                print("Issue in parsing:{}".format(dataset))
                continue

        # Refinement of minima pdb
        # TODO remove explicit call to 0001

        es_refine_pdb = os.path.join(out_dir, compound, dataset,
                                     "exhaustive_search0001", "es_refine.pdb")

        es_refine_mtz = os.path.join(out_dir, compound, dataset,
                                     "exhaustive_search0001", "es_refine.mtz")

        input_mtz = os.path.join(in_dir, compound, dataset, "refine.mtz")

        dataset_dir = os.path.join("/dls/labxchem/data/2018/lb18145-55/"
                                 "processing/analysis/initial_model",
                                 dataset)

        cifs = [f for f in os.listdir(dataset_dir) if f.endswith(".cif")]

        if len(cifs) == 1:
            input_cif = os.path.join(dataset_dir, cifs[0])
        elif len(cifs) == 0:
            raise Warning, "No cif found"
        else:
            raise Warning, "Multiple cif files found"

        if not os.path.exists(es_refine_pdb):
            try:
                os.chdir(os.path.join(out_dir, compound, dataset))
                os.system("giant.quick_refine input.pdb={} "
                          "input.mtz={} input.cif={} "
                          "dir_prefix='exhaustive_search' "
                          "output.out_prefix='es_refine' ".format(es_pdb,
                                                             input_mtz,
                                                             input_cif))
            except IOError:
                print("Skipping crystal")
                continue

        es_minima_plot_folder = os.path.join(out_dir, compound,
                                             dataset,"Plots")
        # Plotting spider plots of minima pdb
        if not os.path.exists(es_minima_plot_folder):
            score_params = score_phil.extract()
            score_params.input.pdb1 = es_pdb
            score_params.input.mtz1 = input_mtz
            score_params.input.pdb2 = es_refine_pdb
            score_params.input.mtz2 = es_refine_mtz
            score_params.output.out_dir = es_minima_plot_folder
            score_model(score_params)

        # Minima occupancy and B factor for histogram/ scatter summary
        try:
            occ, u_iso, _ = get_minimum_fofc(csv_name)
        except IOError:
            print("Issue in parsing:{}".format(dataset))
            continue


        es_occs.append(occ)
        es_b_fac.append(u_iso_to_b_fac(u_iso))

    params.output.out_dir = os.path.join(out_dir,compound)

    occupancy_histogram_with_exhaustive_search(es_occs,
                                               refine_occs,
                                               protein_name="NUDT22",
                                               compound=compound,
                                               params=params)

    occupancy_b_factor_scatter_plot(es_occs=es_occs,
                                    refined_occs=refine_occs,
                                    es_b_fac=es_b_fac,
                                    refine_mean_b_fac=mean_ligand_b_factor,
                                    refine_std_b_fac=std_ligand_b_fac,
                                    protein_name="NUDT22A",
                                    compound=compound,
                                    params=params)





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
