import numpoy as np

from exhaustive.utils.utils_ccp4 import round_step


def b_to_u_iso(b_fac):
    """ Convert isotropic B factor to u iso"""

    u_iso = b_fac/(8 * np.pi ** 2)
    return u_iso


def u_iso_to_b_fac(u_iso):
    """ Convert u_iso to isotropic B factor """

    b_iso = (8 * np.pi ** 2) * u_iso
    return b_iso


def get_minimum_fofc(csv_name, b_fac=None):
    """
    Get minima in fofc, and return minima and where it occurs

    B factor can be supplied to look at the minima across a
    single b factor value

    :param csv_name:
    :param b_fac:
    :return:
    """
    print(os.getcwd())
    # TODO Remove this dual .csv by cleaning up csv name: issue 59
    if csv_name.endswith(".csv"):
        data = np.genfromtxt(csv_name, delimiter=',', skip_header=0)
    else:
        data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',',
                             skip_header=0)

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
        raise ValueError("|Fo-Fc| data is not in correct format")

    if b_fac is not None:
        set_u_iso = b_to_u_iso(b_fac)
        step = np.unique(occ)[1] - np.unique(occ)[0]
        check_u_iso = round_step(set_u_iso, base=step)
        data_array = np.stack((occ, u_iso, fo_fc))

        occ = data_array[0][(data_array[1] >= check_u_iso)
                            & (data_array[1] <= check_u_iso)]

        u_iso = data_array[1][(data_array[1] >= check_u_iso)
                              & (data_array[1] <= check_u_iso)]

        fo_fc = data_array[2][(data_array[1] >= check_u_iso)
                              & (data_array[1] <= check_u_iso)]

    min_index = np.argmin(fo_fc)

    return occ[min_index], u_iso[min_index], fo_fc[min_index]