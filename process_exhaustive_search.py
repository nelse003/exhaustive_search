import csv
import os
import sys

import iotbx
import libtbx.phil
import numpy as np
from iotbx.pdb import hierarchy

from Repeating_exhaustive_search import get_in_refinement_or_better
from select_occupancy_groups import process_refined_pdb_bound_ground_states

#################################################
master_phil = libtbx.phil.parse("""
input{
    database_path = "/dls/labxchem/data/2018/lb18145-55/processing/database/soakDBDataFile.sqlite"
        .type = path
    csv_name = 'u_iso_occupancy_vary_new_atoms'
        .type = str
}
output{
    out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/occupancy_group_with_refinement"
        .type = str
    minima_csv_name = "min_occ_u_iso_NUDT22_with_refinement"
        .type = str
}
options{
    cat = "cat"
        .type = str
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

    return occ[min_index], u_iso[min_index], fo_fc[min_index]

# TODO Remove this second copy, by sorting out circular references
def u_iso_to_b_fac(u_iso):

    b_iso = (8 * np.pi ** 2) * u_iso ** 2
    return b_iso

def write_minima_pdb(input_pdb,output_pdb,csv_name):

    min_occ, min_u_iso, _ = get_minimum_fofc(csv_name)

    bound_states, ground_states = process_refined_pdb_bound_ground_states(input_pdb)
    pdb_inp = iotbx.pdb.input(input_pdb)
    hier = pdb_inp.construct_hierarchy()

    for chain in hier.only_model().chains():
        for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                for atom in atom_group.atoms():

                    for ground_state in ground_states:
                        num_altlocs = ground_state[1]
                        if ground_state[0][atom.i_seq]:
                            atom.occ = (1 - min_occ) / num_altlocs
                            atom.b = u_iso_to_b_fac(min_u_iso)

                    for bound_state in bound_states:
                        num_altlocs = bound_state[1]
                        if bound_state[0][atom.i_seq]:
                            atom.set_occ(min_occ/num_altlocs)
                            atom.set_b(u_iso_to_b_fac(min_u_iso))


    with open(output_pdb,"w") as file:
        file.write(hier.as_pdb_string(crystal_symmetry=hierarchy.input(input_pdb).crystal_symmetry()))


def get_all_minima(params):

    start_xtal_num = 0
    end_xtal_num = 2000
    prefix = "DCP2B-x"
    xtals = []
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    with open(params.output.minima_csv_name+'.csv', 'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        for xtal_name, pdb, mtz in get_in_refinement_or_better(params):

            if xtal_name in xtals:

                print("Getting u_iso and occupancy @ sum(|fo-fc|) minima, for {}".format(xtal_name))
                if pdb and mtz is not None:

                    assert os.path.exists(pdb), 'PDB File does not exist: {}'.format(pdb)
                    assert os.path.exists(mtz), 'MTZ File does not exist: {}'.format(mtz)

                    os.chdir(os.path.join(params.output.out_dir, xtal_name))
                    try:
                        occ, u_iso, _ = get_minimum_fofc(params.input.csv_name)
                    except IOError:
                        print("Skipping crystal {}".format(xtal_name))
                        continue
                    row = [xtal_name, occ, u_iso]
                    writer.writerow(row)
                    sys.stdout.flush()
                    os.chdir("../..")

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

def run(params):

    # print(params.input.database_path)
    # get_all_minima(params)
    print("AAAAAAAAAAAAA")
    input_pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1787/refine.pdb"
    csv_name = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/NUDT7_Copied_atoms/NUDT7A-x1787/u_iso_occupancy_vary_new_atoms"
    write_minima_pdb(input_pdb, csv_name)


if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:])