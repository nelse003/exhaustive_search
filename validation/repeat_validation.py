import os

import pandas as pd
from giant.jiffies.split_conformations import master_phil as split_phil
from giant.jiffies.split_conformations import run as split_conformations

from exhaustive.utils import master_phil
from repeat_validate import repeat_validate


def list_files(directory, extension):
    return [f for f in os.listdir(directory) if f.endswith("." + extension)]


params = master_phil.extract()

out_dir = (
    "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation_FALZA/"
)
# loop_dir= "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/"
loop_dir = "/dls/labxchem/data/2016/lb13385-61/processing/analysis/initial_model"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# compound_dirs = [os.path.join(loop_dir, compound_dir) for compound_dir in os.listdir(loop_dir)
#                  if os.path.isdir(os.path.join(loop_dir, compound_dir))]

datasets = []

xtals = [
    "FALZA-x0079",
    "FALZA-x0085",
    "FALZA-x0172",
    "FALZA-x0177",
    "FALZA-x0271",
    "FALZA-x0309",
    "FALZA-x0402",
    "FALZA-x0438",
]

# for compound_dir in compound_dirs:

xtal_dirs = [
    os.path.join(loop_dir, xtal_dir)
    for xtal_dir in os.listdir(loop_dir)
    if os.path.isdir(os.path.join(loop_dir, xtal_dir))
]

# compound_name = os.path.basename(compound_dir)

for xtal_dir in xtal_dirs:

    xtal_name = os.path.basename(xtal_dir)

    if xtal_name in xtals:

        compounds = list_files(xtal_dir, "cif")
        compound_name = (list_files(xtal_dir, "cif")[0]).split(".")[0]

        refine_pdb = os.path.join(xtal_dir, "refine.pdb")
        refine_mtz = os.path.join(xtal_dir, "refine.mtz")
        xtal_out_dir = os.path.join(out_dir, compound_name, xtal_name)

        if not os.path.exists(os.path.join(out_dir, compound_name)):
            os.mkdir(os.path.join(out_dir, compound_name))

        if not os.path.exists(os.path.join(out_dir, compound_name, xtal_name)):
            os.mkdir(os.path.join(out_dir, compound_name, xtal_name))

        if os.path.exists(refine_pdb) and os.path.exists(refine_mtz):
            datasets.append(
                (
                    xtal_name,
                    compound_name,
                    xtal_dir,
                    refine_pdb,
                    refine_mtz,
                    xtal_out_dir,
                )
            )
        else:
            continue

print(datasets)

# datasets = []
# for compound_dir in compound_dirs:
#
#     xtal_dirs = [os.path.join(compound_dir,xtal_dir) for xtal_dir in os.listdir(compound_dir)
#                  if os.path.isdir(os.path.join(compound_dir, xtal_dir))]
#
#     compound_name = os.path.basename(compound_dir)
#
#     for xtal_dir in xtal_dirs:
#
#         xtal_name = os.path.basename(xtal_dir)
#         refine_pdb = os.path.join(xtal_dir,"refine.pdb")
#         refine_mtz = os.path.join(xtal_dir,"refine.mtz")
#         xtal_out_dir = os.path.join(out_dir, compound_name, xtal_name)
#
#         if not os.path.exists(os.path.join(out_dir, compound_name)):
#             os.mkdir(os.path.join(out_dir, compound_name))
#
#         if not os.path.exists(os.path.join(out_dir, compound_name, xtal_name)):
#             os.mkdir(os.path.join(out_dir, compound_name, xtal_name))
#
#         if os.path.exists(refine_pdb) and os.path.exists(refine_mtz):
#             datasets.append((xtal_name, compound_name, xtal_dir, refine_pdb, refine_mtz, xtal_out_dir))
#         else:
#             continue

# validation based params

params.exhaustive.options.column_type = "FMODEL"
params.exhaustive.options.generate_mtz = False
params.validate.options.use_qsub = False
params.validate.options.step_simulation = 0.05
params.validate.options.overwrite = True
params.exhaustive.options.step = 0.05
params.settings.processes = 14
params.validate.options.set_b = 40.0

for dataset in datasets:

    print(dataset)

    (
        params.input.xtal_name,
        compound_name,
        params.input.in_path,
        params.input.pdb,
        params.input.mtz,
        params.output.out_dir,
    ) = dataset

    params.validate.input.base_mtz = params.input.mtz
    params.output.log_dir = os.path.join(params.output.out_dir, "logs")

    params.validate.input.ground_state_pdb_path = os.path.join(
        params.input.in_path, "refine.output.ground-state.pdb"
    )
    params.validate.input.bound_state_pdb_path = os.path.join(
        params.input.in_path, "refine.output.bound-state.pdb"
    )

    if (
        not os.path.exists(params.validate.input.ground_state_pdb_path)
        or params.validate.options.overwrite
    ):
        split_params = split_phil.extract()
        split_params.input.pdb = [params.input.pdb]
        split_params.output.suffix_prefix = "output"
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

    params.validate.options.repeat_validate_qsub = False

    print(params.output.out_dir)

    if params.validate.options.repeat_validate_qsub:
        modified_phil = master_phil.format(python_object=params)

        with open(
            os.path.join(params.output.out_dir, "params.txt"), "w+"
        ) as param_file:
            param_file.write(modified_phil.as_str())
        with open(
            os.path.join(params.output.out_dir, "run_repeat_validation.py"), "w+"
        ) as python_file:
            python_file.write("import os, sys\n")
            python_file.write("from libtbx.phil import parse\n")
            python_file.write(
                "scriptpath='/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search'\n"
            )
            python_file.write("sys.path.insert(0, os.path.abspath(scriptpath))\n")
            python_file.write("from phil import master_phil\n")
            python_file.write(
                "from exhaustive.validation.repeat_validate import repeat_validate\n"
            )
            # python_file.write('file = open(os.path.join(\'{}\',"params.txt"))\n'.format(params.output.out_dir))
            # python_file.write('params_string = file.read()\n')
            # python_file.write('file.close()\n')
            # python_file.write('print(params_string)\n')
            python_file.write(
                "user_phil=parse(file_name=os.path.join('{}',\"params.txt\"))\n".format(
                    params.output.out_dir
                )
            )
            python_file.write("working_phil = master_phil.fetch(sources=[user_phil])\n")
            python_file.write("params =  working_phil.extract()\n")
            python_file.write("print(params.output.out_dir)\n")
            python_file.write("repeat_validate(params)\n")

        with open(
            os.path.join(params.output.out_dir, "run_repeat_validation.sh"), "w"
        ) as f:
            f.write("#!/bin/bash\n")
            f.write(
                "source /dls/science/groups/i04-1/software/pandda-update/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n"
            )
            f.write(
                "/dls/science/groups/i04-1/software/pandda-update/ccp4/ccp4-7.0/bin/ccp4-python {}".format(
                    os.path.join(params.output.out_dir, "run_repeat_validation.py")
                )
            )

        # This qsub is failing becuase it can't import libtbx.
        # However libtbx should be provided by ccp4-python call
        # non-qsub submission of

        print(
            "qsub {}".format(
                "$CCP4/bin/ccp4-python {}".format(
                    os.path.join(params.output.out_dir, "run_repeat_validation.py")
                )
            )
        )

        # os.system('qsub {}'.format("$CCP4/bin/ccp4-python {}".format(
        #     os.path.join(params.output.out_dir, "run_repeat_validation.py"))))

    else:
        repeat_validate(params)

validation_summary_dfs = []
for dataset in datasets:

    for ext in ["lig_grid", "per_residue", "convex_hull_buffer_0_0"]:

        (
            params.input.xtal_name,
            compound_name,
            params.input.in_path,
            params.input.pdb,
            params.input.mtz,
            params.output.out_dir,
        ) = dataset

        csv_path = os.path.join(params.output.out_dir, ext, "validation_summary.csv")
        try:
            df = pd.read_csv(csv_path)
            df["method"] = ext
            df["compound"] = compound_name
            df = df.rename(index={0: params.input.xtal_name})
            print(df)
            validation_summary_dfs.append(df)
        except IOError:
            print("{}: Not done".format(params.input.xtal_name))
            continue

df = pd.concat(validation_summary_dfs)
print(df)
df.to_csv(os.path.join(out_dir, "validation_summary_all.csv"))
