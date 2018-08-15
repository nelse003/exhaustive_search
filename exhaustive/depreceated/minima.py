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


def get_all_minima(params):
    start_xtal_num = 0
    end_xtal_num = 2000
    prefix = "DCP2B-x"
    xtals = []
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    with open(params.output.minima_csv_name + '.csv', 'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        for xtal_name, pdb, mtz in get_in_refinement_or_better(params):

            if xtal_name in xtals:

                print("Getting u_iso and occupancy @ sum(|fo-fc|) minima, for {}".format(xtal_name))
                if pdb and mtz is not None:

                    try:
                        assert os.path.exists(pdb), 'PDB File does not exist: {}'.format(pdb)
                    except AssertionError:
                        logging.exception('PDB File does not exist: %s', pdb)
                        raise
                    try:
                        assert os.path.exists(mtz), 'MTZ File does not exist: {}'.format(mtz)
                    except AssertionError:
                        logging.exception('MTZ File does not exist: %s', mtz)
                        raise

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

    idx = np.argpartition(mean_fo_fc, 10)
    print(u_iso[idx[:10]], occ[idx[:10]], mean_fo_fc[idx[:10]])
    u_iso_mean = np.mean(u_iso[idx[:10]])
    u_iso_std = np.std(u_iso[idx[:10]])
    occ_mean = np.mean(occ[idx[:10]])
    occ_std = np.std(occ[idx[:10]])
    mean_fofc_mean = np.mean(mean_fo_fc[idx[:10]])
    mean_fofc_std = np.std(mean_fo_fc[idx[:10]])

    print(u_iso_mean, u_iso_std, occ_mean, occ_std, mean_fofc_mean, mean_fofc_std)