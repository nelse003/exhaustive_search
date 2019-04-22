from giant import grid as grid
from giant.maths.geometry import pairwise_dists

from utils.select_atoms import logging

def get_parameter_from_occupancy_groups(occupancy_groups, parameter_str):
    """
    Extracts a parameter from occupancy groups,
    given a str matching that parameter (altloc, chain, resseq...)

    :param occupancy_groups:
    :param parameter_str:
    :return:
    """

    parameters = []

    for occupancy_group in occupancy_groups:
        for group in occupancy_group:
            for residue_altloc in group:
                if residue_altloc.get("model") == '':
                    parameters.append(residue_altloc.get(parameter_str))
                else:
                    raise Warning("Multiple models are present in pdb file. "
                                  "This is not processable with occupancy "
                                  "group selection")
    if not parameters:
        logging.warning("Parameter may not be recognised,"
                        "as output list is empty")
        raise Warning("Parameter may not be recognised,"
                      "as output list is empty")

    return parameters


def within_rmsd_cutoff(atoms1, atoms2, params):
    """
    Given two groups of atoms determine within a given cutoff
    (supplied via params)

    :param atoms1:
    :param atoms2:
    :param params:
    :return:
    """

    for i in range(0, len(pairwise_dists(atoms1.extract_xyz(),
                                         atoms2.extract_xyz()))):
        if pairwise_dists(atoms1.extract_xyz(),
                          atoms2.extract_xyz())[i][i] \
                < params.select.coincident_cutoff:
            continue
        else:
            return False
    return True