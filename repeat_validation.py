import os

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

start_xtal_num = 909
end_xtal_num = 937
in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/repeat_soaks/2018-05-28/NUDT22_from_occ_group_with_refinement/FMOPL000622a_DSPL"
out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/validation_NUDT22/FMOPL000622a_DSPL"
prefix = "NUDT22A-x"


#validation based params

params.exhaustive.options.column_type = "FMODEL"
params.exhaustive.options.generate_mtz = False
params.validate.options.use_qsub = False
params.validate.options.step_simulation = 0.05
params.validate.options.overwrite = False
params.exhaustive.options.step = 0.05
params.settings.processes = 14

# # copy data to new folder

if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    #os.system('cp -a {}/. {}'.format(in_dir,out_dir))

xtals = ["NUDT22A-x0182"]
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = prefix + "{0:0>4}".format(num)
    xtals.append(xtal_name)

for xtal_name in xtals:

    params.input.xtal_name = xtal_name
    params.input.pdb = os.path.join(os.path.join(in_dir, xtal_name, "refine.pdb"))
    params.input.mtz = os.path.join(os.path.join(in_dir, xtal_name, "refine.mtz"))

    params.validate.input.base_mtz = os.path.join(os.path.join(in_dir,
                                                               xtal_name,
                                                               "refine.mtz"))
    print(params.input.pdb)

    if not os.path.exists(params.input.pdb):
        continue
    if not os.path.exists(params.input.mtz):
        continue

    params.output.out_dir = os.path.join(out_dir, xtal_name)
    params.output.log_dir = os.path.join(out_dir, xtal_name, "logs")
    #params.exhaustive.output.csv_name = os.path.join(params.output.out_dir, "exhaustive_search.csv")

    params.input.in_path = os.path.join(os.path.join(in_dir, xtal_name))
    params.output.out_dir = os.path.join(os.path.join(out_dir, xtal_name))
    params.output.log_dir = os.path.join(params.output.out_dir, "logs")
    params.validate.input.ground_state_pdb_path = os.path.join(
        params.input.in_path, "refine.output.ground-state.pdb")
    params.validate.input.bound_state_pdb_path = os.path.join(
        params.input.in_path, "refine.output.bound-state.pdb")
    params.validate.options.set_b = 40.0

    if not os.path.exists(params.validate.input.ground_state_pdb_path):
        split_params = split_phil.extract()
        split_params.input.pdb = [params.input.pdb]
        split_params.output.suffix_prefix = 'output'
        split_conformations(split_params)

    if not os.path.exists(params.output.out_dir):
        os.mkdir(params.output.out_dir)

    if not os.path.exists(params.output.log_dir):
        os.mkdir(params.output.log_dir)

    validate(params)

    # Removal of existing output files for cctbx fmodel to run
    if os.path.exists(params.output.out_dir):

        for item in os.listdir(params.output.out_dir):
            if item.endswith(".mtz"):
                if not item.startswith("refine"):
                    os.remove(os.path.join(params.output.out_dir, item))