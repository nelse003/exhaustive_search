import numpoy as np

def b_to_u_iso(b_fac):
    """ Convert isotropic B factor to u iso"""

    u_iso = b_fac/(8 * np.pi ** 2)
    return u_iso


def u_iso_to_b_fac(u_iso):
    """ Convert u_iso to isotropic B factor """

    b_iso = (8 * np.pi ** 2) * u_iso
    return b_iso