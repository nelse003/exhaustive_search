import os

import numpy as np

# from iotbx.pdb import hierarchy
from plot_exhaustive_search import scatter_plot


def buffer_validation():

    input_pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/refine.pdb"
    input_mtz = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/refine.mtz"
    out_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/buffer_vary/"
    xtal_name = "NUDT7A-x1237"

    for buffer in np.arange(0.0,0.76,0.25):

        out_dir = os.path.join(out_path, "{}_{}".format(xtal_name, buffer))
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        buffer_txt = str(buffer).replace(".","_")
        sh_file = "{}_{}.sh".format(xtal_name,buffer_txt)
        #
        # file = open(os.path.join(out_dir, sh_file),'w')
        #
        # file.write("#!/bin/bash\n")
        # file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
        # file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
        # file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive_search.py "\
        #            "input.pdb={} input.mtz={} options.buffer={} "\
        #            "xtal_name={} output.out_dir={} \n".format(input_pdb, input_mtz, buffer,xtal_name, out_dir))
        #
        # file.close()

        # os.system("qsub -o {} -e {} {}".format(os.path.join(out_dir,"output.txt"),
        #                                      os.path.join(out_dir,"error.txt")
        #                                     ,os.path.join(out_dir, sh_file)))
        scatter_plot(os.path.join(out_dir,"u_iso_occupancy_vary"))

input_pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/refine.pdb"
input_mtz = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/NUDT7A-x1237.mtz"
xtal_name = "NUDT7A-x1237"

os.chdir("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation")

for occupancy in np.arange(0,1.01,0.05):
    # pdb_in = hierarchy.input(input_pdb)
    #
    # for chain in pdb_in.hierarchy.only_model().chains() :
    #     for residue_group in chain.residue_groups() :
    #         for atom_group in residue_group.atom_groups() :
    #             if atom_group.resname == "LIG":
    #                 for atom in atom_group.atoms() :
    #                     atom.set_occ(occupancy)
    #
    occ_pdb = "NUDT7A-x1237-occ_{}.pdb".format(str(occupancy).replace(".","_"))
    # f = open(occ_pdb, "w")
    # f.write(pdb_in.hierarchy.as_pdb_string(
    #     crystal_symmetry=pdb_in.input.crystal_symmetry()))
    # f.close()
    #
    log_file = "simul_from_refine_log_{}".format(str(occupancy).replace(".","_"))
    mtz_out = "simul_occ_{}.mtz".format(str(occupancy).replace(".","_"))

    os.system("ccp4-python ../simulate_experimental_data.py input.xray_data.file_name={} "
              "model.file_name={} input.xray_data.label=\"I(+),SIGI(+),I(-),SIGI(-),merged\" "
              "output.logfile={} output.hklout={}".format(input_mtz, occ_pdb, log_file,mtz_out))

    exit()
    # input_mtz = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/simul_{}.mtz".format(str(occupancy).replace(".","_"))
    #
    # out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/"
    # csv_name = "occ_u_iso_{}".format(str(occupancy).replace(".","_"))
    #
    # sh_file = "{}_{}.sh".format(xtal_name, occupancy)
    #
    # file = open(os.path.join(out_dir, sh_file),'w')
    #
    # file.write("#!/bin/bash\n")
    # file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
    # file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
    #
    # file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive_search.py"
    #            " input.pdb={} input.mtz={} output.out_dir={} xtal_name={} "
    #            "options.csv_name={}".format(input_pdb,input_mtz,out_dir,xtal_name,csv_name))
    # file.close()
    #
    # os.system("qsub -o {} -e {} {}".format(os.path.join(out_dir,"output_{}.txt".format(str(occupancy).replace(".","_"))),
    #                                        os.path.join(out_dir,"error_{}.txt".format(str(occupancy).replace(".","_"))),
    #                                        os.path.join(out_dir, sh_file)))

# args = [input_pdb, input_mtz]
#
# pdb_inp = iotbx.pdb.input(input_pdb)
#
# inputs = mmtbx.utils.process_command_line_args(args = args)
# ph = pdb_inp.construct_hierarchy()
# xrs = ph.extract_xray_structure(
#     crystal_symmetry = inputs.crystal_symmetry)
# reflection_files = inputs.reflection_files
# rfs = reflection_file_utils.reflection_file_server(
#     crystal_symmetry=inputs.crystal_symmetry,
#     force_symmetry=True,
#     reflection_files=reflection_files,
#     err=StringIO())
# determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
#     reflection_file_server=rfs,
#     keep_going=True,
#     log=StringIO())
# r_free_flags = determined_data_and_flags.r_free_flags
# f_obs = determined_data_and_flags.f_obs
# crystal_gridding = f_obs.crystal_gridding(
#     d_min             = f_obs.d_min(),
#     symmetry_flags    = maptbx.use_space_group_symmetry,
#     resolution_factor = 1./4)
#
# mask_params = mmtbx.masks.mask_master_params.extract()
# mask_params.ignore_hydrogens=False
# mask_params.ignore_zero_occupancy_atoms=False
# fmodel = mmtbx.f_model.manager(
#     r_free_flags   = r_free_flags
#     mask_params    = mask_params,
#     xray_structure = xrs)
# fmodel.update_all_scales()
#
# fcfc_map, fcfc  = compute_maps(fmodel=fmodel,crystal_gridding=crystal_gridding,map_type="mFo-DFc")
# mtz_dataset = fcfc.as_mtz_dataset(column_root_label="FOFCWT")
# mtz_object = mtz_dataset.mtz_object()
# mtz_object.write(file_name="testing_{}_{}.mtz".format(bound_occupancy, u_iso))