import csv

import libtbx.phil
import numpy as np

from Repeating_exhaustive_search import get_in_refinement_or_better

#################################################
master_phil = libtbx.phil.parse("""
input{
    database_path = "/dls/labxchem/data/2018/lb18145-55/processing/database/soakDBDataFile.sqlite"
        .type = path
    csv_name = 'u_iso_occupancy_vary_new_atoms'
        .type = str
}
output{
    out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/occupancy_group_test"
        .type = str
    minima_csv_name = "min_occ_u_iso_all"
        .type = str
}
options{

}
""", process_includes=True)
###################################################

def get_minimum_fofc(csv_name):
    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)

    # If four column data from multiple ligand
    if len(data[0]) == 4:
        occ = data[:, 0]
        u_iso = data[:, 2]
        fo_fc = data[:, 3]
    elif len(data[0]) == 3:
        occ = data[:, 0]
        u_iso = data[:, 1]
        fo_fc = data[:, 2]
    else:
        print("Data is not in correct format")
    # b_iso = (8 * np.pi ** 2) * u_iso ** 2

    # if three column data

    min_index = np.argmin(fo_fc)
    return occ[min_index], u_iso[min_index]


def get_all_minima(params):
    with open(params.output.minima_csv_name, 'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        for xtal_name, pdb, mtz in get_in_refinement_or_better(params):
            print("Getting u_iso and occupancy @ sum(|fo-fc|) minima, for {}".format(xtal_name))
            if pdb and mtz is not None:
                try:
                    assert os.path.exists(pdb), 'PDB File does not exist: {}'.format(pdb)
                    assert os.path.exists(mtz), 'MTZ File does not exist: {}'.format(mtz)
                    os.chdir(os.path.join(params.output.out_dir, xtal_name))
                    occ, u_iso = get_minimum_fofc(params.input.csv_name)
                    row = [xtal_name, occ, u_iso]
                    writer.writerow(row)
                    sys.stdout.flush()
                    os.chdir("../..")
                except:
                    print("Minima processing failed on xtal: {}".format(xtal_name))
                    continue
            else:
                print("Path to PDB & MTZ file is likely incorrect")

def minima_flatness(csv_path):
    """ Estiamte the flatness of a minima"""
    data = np.genfromtxt('{}.csv'.format(csv_path), delimiter=',', skip_header=0)
    occ = data[:, 0]
    u_iso = data[:, 2]
    mean_fo_fc = data[:, 3]

    idx = np.argpartition(mean_fo_fc,10)
    print(u_iso[idx[:10]],occ[idx[:10]],mean_fo_fc[idx[:10]])
    u_iso_mean = np.mean(u_iso[idx[:10]])
    u_iso_std = np.std(u_iso[idx[:10]])
    occ_mean = np.mean(occ[idx[:10]])
    occ_std = np.std(occ[idx[:10]])
    mean_fofc_mean = np.mean(mean_fo_fc[idx[:10]])
    mean_fofc_std = np.std(mean_fo_fc[idx[:10]])

    print(u_iso_mean,u_iso_std, occ_mean, occ_std, mean_fofc_mean, mean_fofc_std)

def check_whether_ground_and_bound_states_exist():
    """ Check whetehr both ground and bound states exist"""
    pass

# minima_flatness('NUDT22A-x1058/NUDT22A-x1058/u_iso_occupancy_vary')
# minima_flatness('occupancy_group_test/NUDT22A-x1058/u_iso_occupancy_vary')

get_all_minima(params = master_phil.extract())