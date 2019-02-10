from exhaustive.exhaustive.utils.utils import read_occupancies_from_refmac_log
from exhaustive.exhaustive.plotting.plot import plot_occupancy_convergence

occ_conv_df = read_occupancies_from_refmac_log("/home/enelson/Downloads/output.quick-refine.log")
plot_occupancy_convergence(occ_conv_df=occ_conv_df,
                           plot_filename="/home/enelson/Desktop/occ_conv.png")
