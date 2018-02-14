import numpy as np

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
