from exhaustive.utils.phil import master_phil
from exhaustive.utils.phil import prepare_validate_phil
from exhaustive.plotting.plot import plot_2d_occ_b_validation

params = master_phil.extract()

plot_2d_occ_b_validation(
    start_occ=0.05,
    end_occ=1,
    step=0.1,
    set_b=40.0,
    dataset_prefix='FALZA-x0085',
    out_dir="/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/test/output/per_res",
    params=params,
)
