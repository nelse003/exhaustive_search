from giant import grid as grid
from giant.maths.geometry import pairwise_dists

from utils.select_atoms import logging


def get_occupancy_group_grid_points(pdb, bound_states, ground_states,
                                    params):
    """Produce cartesian points related to occupancy groups.

    Get cartesian points that correspond to atoms involved in the
    occupancy groups (as in multi-state.restraints.params)

    :param pdb: Input PDB file
    :type
    :param params: Working phil parameters
    :type
    :return: occupancy_group_cart_points: The cartesian points involved
    in the bound and ground states as a list
    """

    logging.info("For all bound and ground states, "
                 "select cartesian grid points for each altloc/residue \n"
                 "involved in occupancy groups. A buffer of {} Angstrom \n"
                 "is applied to minimal and maximal grid points,"
                 "with a grid seperation of {}.\n".format(
        params.exhaustive.options.buffer,
        params.exhaustive.options.grid_spacing))

    states = bound_states + ground_states

    pdb_in = iotbx.pdb.hierarchy.input(pdb)
    pdb_atoms = pdb_in.hierarchy.atoms()

    occupancy_group_cart_points = flex.vec3_double()
    for state in states:
        selection = state[0]
        selected_atoms = pdb_atoms.select(selection)
        sites_cart = selected_atoms.extract_xyz()
        grid_min = flex.double([s - params.exhaustive.options.buffer
                                for s in sites_cart.min()])
        grid_max = flex.double([s + params.exhaustive.options.buffer
                                for s in sites_cart.max()])
        grid_from_selection = grid.Grid(
            grid_spacing=params.exhaustive.options.grid_spacing,
            origin=tuple(grid_min),
            approx_max=tuple(grid_max))

        logging.debug(grid_from_selection.summary())

        occupancy_group_cart_points = occupancy_group_cart_points.concatenate(
            grid_from_selection.cart_points())

    logging.info("Number of cartesian points to calculate "
                 "|Fo-Fc| over: {}".format(len(occupancy_group_cart_points)))

    return occupancy_group_cart_points


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