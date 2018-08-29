import os
import numpy as np

from exhaustive.validation.validation import run as validate
from phil import master_phil
from giant.jiffies.split_conformations import master_phil as split_phil
from giant.jiffies.split_conformations import run as split_conformations

params =  master_phil.extract()

# start_xtal_num = 1905
# end_xtal_num = 1915
# in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios"
# out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation_covalent_ratios"
# prefix = "NUDT7A-x"

# start_xtal_num = 909
# end_xtal_num = 937
# in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/FMOPL000622a_DSPL"
# out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation_NUDT22/FMOPL000622a_DSPL"
# prefix = "NUDT22A-x"


out_dir =  "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation_NUDT22/"
loop_dir= "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/"

compound_dirs = [os.path.join(loop_dir, compound_dir) for compound_dir in os.listdir(loop_dir)
                 if os.path.isdir(os.path.join(loop_dir, compound_dir))]

datasets = []
for compound_dir in compound_dirs:

    xtal_dirs = [os.path.join(compound_dir,xtal_dir) for xtal_dir in os.listdir(compound_dir)
                 if os.path.isdir(os.path.join(compound_dir, xtal_dir))]

    compound_name = os.path.basename(compound_dir)

    for xtal_dir in xtal_dirs:

        xtal_name = os.path.basename(xtal_dir)
        refine_pdb = os.path.join(xtal_dir,"refine.pdb")
        refine_mtz = os.path.join(xtal_dir,"refine.mtz")
        xtal_out_dir = os.path.join(out_dir, compound_name, xtal_dir)

        if os.path.exists(refine_pdb) and os.path.exists(refine_mtz):
            datasets.append((xtal_name, xtal_dir, refine_pdb, refine_mtz, xtal_out_dir))
        else:
            continue


#validation based params

params.exhaustive.options.column_type = "FMODEL"
params.exhaustive.options.generate_mtz = False
params.validate.options.use_qsub = False
params.validate.options.step_simulation = 0.05
params.validate.options.overwrite = True
params.exhaustive.options.step = 0.05
params.settings.processes = 20

# # # copy data to new folder
#
# if not os.path.exists("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation_NUDT22/"):
#     os.mkdir("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation_NUDT22/")
#
# if not os.path.exists(out_dir):
#     os.mkdir(out_dir)
#     #os.system('cp -a {}/. {}'.format(in_dir,out_dir))
#
# xtals = ["NUDT22A-x0182"]
# for num in range(start_xtal_num, end_xtal_num + 1):
#     xtal_name = prefix + "{0:0>4}".format(num)
#     xtals.append(xtal_name)
#
# for xtal_name in xtals:

for dataset in datasets:

    print(dataset)

    (params.input.xtal_name, params.input.in_path,  params.input.pdb, params.input.mtz, params.output.out_dir) = dataset
    params.validate.input.base_mtz = params.validate.input.base_mtz
    params.output.log_dir = os.path.join(params.output.out_dir, "logs")

    # params.input.xtal_name = xtal_name
    # params.input.pdb = os.path.join(os.path.join(in_dir, xtal_name, "refine.pdb"))
    # params.input.mtz = os.path.join(os.path.join(in_dir, xtal_name, "refine.mtz"))
    #
    # params.validate.input.base_mtz = os.path.join(os.path.join(in_dir,
    #                                                            xtal_name,
    #                                                            "refine.mtz"))
    # print(params.input.pdb)
    #
    # if not os.path.exists(params.input.pdb):
    #     continue
    # if not os.path.exists(params.input.mtz):
    #     continue

    # params.output.out_dir = os.path.join(out_dir, xtal_name)
    # params.output.log_dir = os.path.join(out_dir, xtal_name, "logs")
    #params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")

    params.validate.input.ground_state_pdb_path = os.path.join(
        params.input.in_path, "refine.output.ground-state.pdb")
    params.validate.input.bound_state_pdb_path = os.path.join(
        params.input.in_path, "refine.output.bound-state.pdb")

    if not os.path.exists(params.validate.input.ground_state_pdb_path) or params.validate.options.overwrite:

        split_params = split_phil.extract()
        split_params.input.pdb = [params.input.pdb]
        split_params.output.suffix_prefix = 'output'
        split_params.options.reset_occupancies = True
        split_conformations(split_params)

    if not os.path.exists(params.output.out_dir):
        os.mkdir(params.output.out_dir)

    if not os.path.exists(params.output.log_dir):
        os.mkdir(params.output.log_dir)


    # Removal of existing output files for cctbx fmodel to run
    if os.path.exists(params.output.out_dir):

        for item in os.listdir(params.output.out_dir):
            if item.endswith(".mtz"):
                if not item.startswith("refine"):
                    os.remove(os.path.join(params.output.out_dir, item))

    # Loop over set B. Fails at setB =100

    # for set_b in np.arange(20,120,5):
    #     params.validate.options.set_b = set_b
    #     params.output.out_dir = os.path.join(os.path.join(out_dir, xtal_name, "set_b_{}".format(str(set_b).replace(".","_"))))
    #
    #     if os.path.exists(params.output.out_dir):
    #
    #         for item in os.listdir(params.output.out_dir):
    #             if item.endswith(".mtz"):
    #                 if not item.startswith("refine"):
    #                     os.remove(os.path.join(params.output.out_dir, item))
    #
    #     validate(params)

    params.validate.options.set_b = 40

    # params.output.out_dir = os.path.join(os.path.join(out_dir, xtal_name, "per_residue"))
    # if os.path.exists(params.output.out_dir):
    #     for item in os.listdir(params.output.out_dir):
    #         if item.endswith(".mtz"):
    #             if not item.startswith("refine"):
    #                 os.remove(os.path.join(params.output.out_dir, item))
    #params.exhaustive.options.per_residue = True

    validate(params)
    exit()

    # Atom selection isnt working? Something to do with merging in validate?

    # params.exhaustive.options.convex_hull = False
    # atom_points_sel_string = "(chain B and altid C and resid 1) or (chain B and altid D resid 1)"
    # params.output.out_dir = os.path.join(
    #     os.path.join(out_dir, xtal_name, "lig_atoms"))
    # params.exhaustive.options.ligand_atom_points = True
    # validate(params)
    #
    # params.exhaustive.options.ligand_atom_points = False
    # params.exhaustive.options.ligand_grid_points = True
    # params.output.out_dir = os.path.join(
    #     os.path.join(out_dir, xtal_name, "lig_grid"))
    # validate(params)

    params.validate.options.set_b = 40

    for buffer in np.arange(0, 3, 0.5):
        params.exhaustive.options.convex_hull_ignore_nearest = False
        params.output.out_dir = os.path.join(os.path.join(out_dir, xtal_name, "test_convex_hull_buffer_{}".format(str(buffer).replace(".","_"))))
        params.exhaustive.options.convex_hull=True
        params.exhaustive.options.buffer=buffer

        if os.path.exists(params.output.out_dir):

            for item in os.listdir(params.output.out_dir):
                if item.endswith(".mtz"):
                    if not item.startswith("refine"):
                        os.remove(os.path.join(params.output.out_dir, item))

        validate(params)

        params.exhaustive.options.convex_hull_ignore_nearest =True
        params.output.out_dir = os.path.join(
            os.path.join(out_dir, xtal_name, "convex_hull_ignore_nearest_buffer_{}".format(str(buffer).replace(".", "_"))))

        if os.path.exists(params.output.out_dir):

            for item in os.listdir(params.output.out_dir):
                if item.endswith(".mtz"):
                    if not item.startswith("refine"):
                        os.remove(os.path.join(params.output.out_dir, item))

        validate(params)

        params.output.out_dir = os.path.join(
            os.path.join(out_dir, xtal_name, "buffer_{}".format(str(buffer).replace(".", "_"))))
        params.exhaustive.options.convex_hull = False

        if os.path.exists(params.output.out_dir):

            for item in os.listdir(params.output.out_dir):
                if item.endswith(".mtz"):
                    if not item.startswith("refine"):
                        os.remove(os.path.join(params.output.out_dir, item))

        validate(params)


