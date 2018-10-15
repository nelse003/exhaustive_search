from __future__ import division, print_function

import datetime
import os
import shutil

import pandas as pd
from giant.jiffies.score_model import master_phil as score_phil
from giant.jiffies.score_model import run as score_model
from giant.jiffies.split_conformations import master_phil as split_phil
from giant.jiffies.split_conformations import run as split_conformations

from exhaustive.exhaustive.exhaustive import run as exhaustive
from exhaustive.exhaustive.plotting.plot import occupancy_b_factor_scatter_plot
from exhaustive.exhaustive.plotting.plot import occupancy_histogram_with_exhaustive_search
from exhaustive.exhaustive.plotting.plot import plot_edstats_across_soaks
from exhaustive.exhaustive.utils.utils import collate_edstats_scores
from exhaustive.exhaustive.utils.utils import get_minimum_fofc
from exhaustive.exhaustive.utils.utils import get_occ_b
from exhaustive.exhaustive.utils.utils import u_iso_to_b_fac
from exhaustive.exhaustive.utils.utils import datasets_from_compound
#from exhaustive.process.minima import write_minima_pdb
from phil import master_phil
from phil import check_input_files

def get_cif_file_from_dataset(dataset_dir, preferred_cif=None):

    cifs = [f for f in os.listdir(dataset_dir) if f.endswith(".cif")]

    if preferred_cif:
        if os.path.exists(preferred_cif):
            return preferred_cif
        else:
            raise Warning, "Prefered Cif file doesn't exist"
    if len(cifs) == 1:
        input_cif = os.path.join(dataset_dir, cifs[0])
        return input_cif
    elif len(cifs) == 0:
        raise Warning, "No cif found"
    else:
        raise Warning, "Multiple cif files found: {}".format(dataset_dir)



# TODO Add logging statements
def run_es_many_xtals(xtals, in_dir, out_dir, params):

    rejects = []

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for xtal_name in xtals:
        params.input.xtal_name = xtal_name
        params.input.in_path = os.path.join(in_dir,xtal_name)
        params.input.pdb = os.path.join(params.input.in_path, "refine.pdb")
        print(params.input.pdb)
        params.output.out_dir = os.path.join(out_dir, xtal_name)

        if not os.path.exists(params.output.out_dir):
            os.mkdir(params.output.out_dir)

        params.output.log_dir = os.path.join(params.output.out_dir, "logs")

        if not os.path.exists(params.output.log_dir):
            os.mkdir(params.output.log_dir)

        # TODO Standardise the output name?
        params.exhaustive.output.csv_name = params.input.xtal_name \
                                            + "_exhaustive_search_occ_u_iso.csv"

        if os.path.exists(os.path.join(params.output.out_dir,
                          params.exhaustive.output.csv_name)):
            continue

        params.input.mtz = os.path.join(params.input.in_path,"refine.mtz")

        try:
            check_input_files(params)
        except:
            rejects.append(xtal_name)
            continue
        try:
            exhaustive(params)
        except UnboundLocalError:
            rejects.append(xtal_name)
            continue

    print(rejects)

def process_exhaustive_search(compound_codes,
                              initial_model_dir,
                              in_dir,
                              out_dir,
                              protein_name,
                              preferred_cif=None):

    protein_prefix = protein_name + "-x"

    for compound in compound_codes:

        es_occ_b = []

        for dataset in datasets_from_compound(protein_prefix,
                                              compound_folder=os.path.join(
                                                  in_dir,compound)):
            # Define paths
            es_csv = os.path.join(out_dir, compound, dataset,
                                  dataset + "_exhaustive_search_occ_u_iso.csv")
            refine_pdb = os.path.join(in_dir, compound, dataset, "refine.pdb")
            es_pdb = os.path.join(out_dir, compound, dataset, "es_minima.pdb")
            es_refine_pdb = os.path.join(out_dir, compound, dataset,
                                         "exhaustive_search0001",
                                         "es_refine.pdb")
            es_refine_mtz = os.path.join(out_dir, compound, dataset,
                                         "exhaustive_search0001",
                                         "es_refine.mtz")
            input_mtz = os.path.join(in_dir, compound, dataset, "refine.mtz")
            # input_mtz = os.path.join(initial_model_dir, dataset,
            #                          dataset + ".free.mtz")
            input_mtz = os.path.join(initial_model_dir, dataset,
                                     "refine.mtz")
            input_pdb = os.path.join(initial_model_dir, dataset,
                                     "refine.pdb")
            print(input_mtz)
            print(input_pdb)
            print("---------------------")

            dataset_dir = os.path.join(initial_model_dir, dataset)

            # Generate split conformations
            if not os.path.exists(os.path.join(in_dir, compound, dataset,
                                               "refine.split.ground-state.pdb")):
                split_params = split_phil.extract()
                split_params.input.pdb = [
                    os.path.join(in_dir, compound, dataset,
                                 "refine.pdb")]
                split_params.output.log = "cat.log"
                print(split_params.input.pdb)
                print(split_phil.format(python_object=split_params).as_str())
                try:
                    split_conformations(split_params)
                except IOError:
                    print("Split confs: Issue in parsing:{}".format(dataset))
                    continue

            # Write Minima pdb
            if not os.path.exists(es_pdb):
                try:
                    write_minima_pdb(input_pdb=refine_pdb,
                                     output_pdb=es_pdb,
                                     csv_name=es_csv,
                                     params=params)
                except IOError:
                    print("Issue in parsing:{}".format(dataset))
                    continue

            input_cif = get_cif_file_from_dataset(dataset_dir, preferred_cif)

            # Refinement of minima pdb
            # TODO remove explicit call to 0001
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
                                                 dataset, "Plots")
            # Plotting spider plots of minima pdb
            if not os.path.exists(es_minima_plot_folder):
                score_params = score_phil.extract()
                score_params.input.pdb1 = input_pdb
                score_params.input.mtz1 = input_mtz
                score_params.input.pdb2 = es_refine_pdb
                score_params.input.mtz2 = es_refine_mtz
                score_params.output.out_dir = es_minima_plot_folder
                try:
                    score_model(score_params)
                except:
                    print("skipping edstats for {}".format(dataset))

            # Minima occupancy and B factor for histogram/ scatter summary
            try:
                occ, u_iso, _ = get_minimum_fofc(es_csv)
            except IOError:
                print("Issue in parsing:{}".format(dataset))
                continue

            es_occ_b.append([dataset,occ,u_iso_to_b_fac(u_iso)])


        print(es_occ_b)
        es_occ_df = pd.DataFrame(data=es_occ_b, columns=['dataset',
                                                        'es_occupancy',
                                                        'es_b_fac'])

        refine_occ_df = get_occ_b(
            refinement_dir=os.path.join(in_dir, compound),
            lig_chain="B",
            pdb_name="refine.split.bound-state.pdb")

        occ_df = pd.merge(es_occ_df,refine_occ_df, on='dataset', how='outer')

        params.output.out_dir = os.path.join(out_dir, compound)

        occupancy_histogram_with_exhaustive_search(occ_df,
                                                   protein_name=protein_name,
                                                   compound=compound,
                                                   params=params)
        occupancy_b_factor_scatter_plot(occ_df,
                                        protein_name=protein_name,
                                        compound=compound,
                                        params=params)

        edstats_df = collate_edstats_scores(protein_prefix=protein_prefix,
                                            compound_folder=os.path.join(
                                                out_dir, compound))
        if edstats_df is not None:
            plot_edstats_across_soaks(edstats_df=edstats_df,
                                      compound_folder=os.path.join(out_dir,
                                                                   compound),
                                      compound=compound,
                                      protein_name=protein_name,
                                      title_suffix="Exhaustive minima")


############
# NUDT7A covalent
############

params =  master_phil.extract()

covalent_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios"

params.output.out_dir = covalent_dir

xtal_dirs = [os.path.join(covalent_dir, xtal_dir)
             for xtal_dir in os.listdir(covalent_dir)
             if os.path.isdir(os.path.join(covalent_dir, xtal_dir))
             and not xtal_dir.endswith("LIG_CYS")]

for xtal_dir in xtal_dirs:
    if not os.path.exists(os.path.join("refine.split.bound-state.pdb")):
        os.chdir(xtal_dir)
        os.system("giant.split_conformations refine.pdb")

refine_occ_df = get_occ_b(
    refinement_dir=covalent_dir,
    lig_chain="E",
    pdb_name="refine.split.bound-state.pdb")

exhaustive_search_csvs = [os.path.join(covalent_dir,
                                       xtal_dir,
                                       xtal_dir.rstrip('_LIG_CYS'),
                                       "exhaustive_search.csv")
                          for xtal_dir in os.listdir(covalent_dir)
                          if os.path.isdir(os.path.join(covalent_dir,
                                                        xtal_dir))
                          and xtal_dir.endswith("LIG_CYS")]


es_occ_b = []
for csv in exhaustive_search_csvs:
    if os.path.exists(csv):
        occ, u_iso, _ = get_minimum_fofc(csv)
        es_occ_b.append([os.path.basename(os.path.dirname(csv)).rstrip('_LIG_CYS'),
                         occ,
                         u_iso_to_b_fac(u_iso)])

es_occ_df = pd.DataFrame(data=es_occ_b, columns=['dataset',
                                                'es_occupancy',
                                                'es_b_fac'])

occ_df = pd.merge(es_occ_df, refine_occ_df, on='dataset', how='outer')


occupancy_histogram_with_exhaustive_search(occ_df,
                                           protein_name="NUDT7A",
                                           compound="NUOOA000181a",
                                           params=params)

occupancy_b_factor_scatter_plot(occ_df,
                                protein_name="NUDT7A",
                                compound="NUOOA000181a",
                                params=params)

exit()
#########################


# Setting common parameters
common_out_dir = os.path.join("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/",
                        datetime.datetime.now().strftime('%Y-%m-%d')) #_%H-%M-%S

if not os.path.exists(common_out_dir):
    os.mkdir(common_out_dir)

params.exhaustive.options.generate_mtz = False
params.exhaustive.options.step = 0.01
params.settings.processes = 36

######################################################
# NUDT22 parameters
######################################################

in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/occupancy_group_with_refinement"
out_dir = os.path.join(common_out_dir, "NUDT22_from_occ_group_with_refinement")
initial_model_dir = "/dls/labxchem/data/2018/lb18145-55/processing/"\
                    "analysis/initial_model"

#####################################################
# NUDT22 Run Exhaustive search on all xtals
#####################################################

start_xtal_num = 909
end_xtal_num = 1058
protein_name =  "NUDT22A"
protein_prefix = protein_name + "-x"
NUDT22_xtals = ['NUDT22A-x0182','NUDT22A-x0243', 'NUDT22A-x0421','NUDT22A-x0391']
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = protein_prefix + "{0:0>4}".format(num)
    NUDT22_xtals.append(xtal_name)

out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/" \
          "NUDT22A"


#TODO Check local folder

# run_es_many_xtals(NUDT22_xtals, in_dir, out_dir, params)

######################################################
# NUDT22 Separate crystals into compound sets
######################################################
xtals_with_compound = {}
# # NUDT22A FMOPL00622a
# from DSI poised library (x0938-x0976)
compound_code = 'FMOPL000622a_DSPL'
xtals_with_compound['NUDT22A-x0182'] =  compound_code
start_xtal_num = 909
end_xtal_num = 937
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = protein_prefix + "{0:0>4}".format(num)
    xtals_with_compound[xtal_name] = compound_code

# # NUDT22A FMOPL00622a
# from DSI poised library (x0938-x0976)
compound_code = 'FMOPL000622a_DSI_poised'
start_xtal_num = 938
end_xtal_num = 976
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = protein_prefix + "{0:0>4}".format(num)
    xtals_with_compound[xtal_name] = compound_code

# NUDT22A 133725a x0421, x1040 - x1059

start_xtal_num = 1040
end_xtal_num = 1059
compound_code = '133725a'
xtals_with_compound['NUDT22A-x0421'] = compound_code
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = protein_prefix + "{0:0>4}".format(num)
    xtals_with_compound[xtal_name] = compound_code

# NUDT22A 13369a x0243, x977 to x1008

start_xtal_num = 977
end_xtal_num = 1008
compound_code = '13369a'
xtals_with_compound['NUDT22A-x0243'] = compound_code
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = protein_prefix + "{0:0>4}".format(num)
    xtals_with_compound[xtal_name] = compound_code

# NUDT22A 13663a x0391, x1009 to x1039

compound_code = '13363a'
xtals_with_compound['NUDT22A-x0391']= compound_code
start_xtal_num = 1009
end_xtal_num = 1039
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = protein_prefix + "{0:0>4}".format(num)
    xtals_with_compound[xtal_name] = compound_code

compound_codes=set()
for compounds in xtals_with_compound.values():
    compound_codes.add(compounds)

for compound in compound_codes:
    if not os.path.exists(os.path.join(in_dir,compound)):
        os.mkdir(os.path.join(in_dir,compound))

for xtal, compound_code in xtals_with_compound.iteritems():
    if os.path.exists(os.path.join(in_dir,xtal)):
        shutil.move(os.path.join(in_dir,xtal),os.path.join(in_dir,
                                                           compound_code))

process_exhaustive_search(compound_codes, initial_model_dir, in_dir, out_dir,
                          protein_name)

#######################################################
# NUDT7 copied atoms (OX210)
# Run exhaustive search
#######################################################

NUDT7_copied_dir = ("/dls/science/groups/i04-1/elliot-dev/Work/"
                   "exhaustive_search_data/NUDT7_Copied_atoms")
NUDT7_copied_xtals = []

for dataset in datasets_from_compound(protein_prefix,
                                      compound_folder=NUDT7_copied_dir):

    NUDT7_copied_xtals.append(NUDT7_dataset)

# run_es_many_xtals(xtals=NUDT7_copied_xtals,
#                   in_dir=NUDT7_copied_dir,
#                   out_dir="/dls/science/groups/i04-1/elliot-dev/Work/"
#                           "exhaustive_search_data/repeat_soaks/2018-05-27/"
#                           "NUDT7_copied",
#                   params=params)

################################################################
# NUDT7 copied atoms (OX210)
#
# Plot the copied atoms refinement and exhaustive search occupancy
################################################################

NUDT7_compound_codes = ["OX210"]
NUDT7_initial_model_dir = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model"

NUDT7_in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/"\
          "exhaustive_search_data/NUDT7_Copied_atoms/"

NUDT7_out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/"\
          "exhaustive_search_data/repeat_soaks/2018-05-27/NUDT7_copied"

process_exhaustive_search(compound_codes=NUDT7_compound_codes,
                          initial_model_dir=NUDT7_initial_model_dir,
                          in_dir=NUDT7_in_dir,
                          out_dir=NUDT7_out_dir,
                          protein_name="NUDT7A",
                          preferred_cif=
"/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"\
"NUDT7_Copied_atoms/OX210/NUDT7A-x1787/OX-210.cif")

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

