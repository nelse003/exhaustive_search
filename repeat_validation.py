import os
import numpy as np

from exhaustive.validation.validation import run as validate
from phil import master_phil
from giant.jiffies.split_conformations import master_phil as split_phil
from giant.jiffies.split_conformations import run as split_conformations
from exhaustive.exhaustive.utils.convex_hull import atom_points_from_sel_string
from exhaustive.exhaustive.utils.utils import process_validation_csvs

def repeat_validate(params):

    print(params.output.out_dir)

    # Issue: Can't do selection of lig atoms for each case, can we do it by LIG naming

    # Ligand grid (by convex hull of ligand atoms)
    params.exhaustive.options.ligand_grid_points = True
    params.output.out_dir = os.path.join(params.output.out_dir, "lig_grid")
    validate(params)
    # Add plotting of residue selection

    # Reset
    params.exhaustive.options.ligand_grid_points = False

    # Per residue selection of atoms
    params.exhaustive.options.per_residue = True
    params.output.out_dir = os.path.join(params.output.out_dir, "per_residue")
    if os.path.exists(params.output.out_dir):
        for item in os.listdir(params.output.out_dir):
            if item.endswith(".mtz"):
                if not item.startswith("refine"):
                    os.remove(os.path.join(params.output.out_dir, item))
    validate(params)
    # Add plotting of residue selection


    # Buffer range (convex hull around occupancy group)
    for buffer in np.arange(0, 2, 0.5):
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

params =  master_phil.extract()

out_dir =  "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation_NUDT22/"
loop_dir= "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/"

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

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
        xtal_out_dir = os.path.join(out_dir, compound_name, xtal_name)

        if not os.path.exists(os.path.join(out_dir, compound_name)):
            os.mkdir(os.path.join(out_dir, compound_name))

        if not os.path.exists(os.path.join(out_dir, compound_name, xtal_name)):
            os.mkdir(os.path.join(out_dir, compound_name, xtal_name))

        if os.path.exists(refine_pdb) and os.path.exists(refine_mtz):
            datasets.append((xtal_name, xtal_dir, refine_pdb, refine_mtz, xtal_out_dir))
        else:
            continue


#validation based params

params.exhaustive.options.column_type = "FMODEL"
params.exhaustive.options.generate_mtz = False
params.validate.options.use_qsub = False
params.validate.options.step_simulation = 0.05
params.validate.options.overwrite = False
params.exhaustive.options.step = 0.05
params.settings.processes = 20
params.validate.options.set_b = 40.0

for dataset in datasets:

    print(dataset)

    (params.input.xtal_name, params.input.in_path,  params.input.pdb, params.input.mtz, params.output.out_dir) = dataset
    params.validate.input.base_mtz = params.input.mtz
    params.output.log_dir = os.path.join(params.output.out_dir, "logs")

    params.validate.input.ground_state_pdb_path = os.path.join(
        params.input.in_path, "refine.output.ground-state.pdb")
    params.validate.input.bound_state_pdb_path = os.path.join(
        params.input.in_path, "refine.output.bound-state.pdb")

    print(atom_points_from_sel_string(params, "resseq LIG"))
    exit()

    # Turn into function, move to after repeat_validate

    # min_fofcs, min_occs, min_b_facs, fofcs, occs, b_facs = \
    #     process_validation_csvs(params.validate.options.start_simul_occ,
    #                  params.validate.options.end_simul_occ,
    #                  step=params.validate.options.step_simulation,
    #                  set_b=params.validate.options.set_b,
    #                  out_dir=params.output.out_dir,
    #                  params=params)

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


    repeat_validate(params)

    exit()






