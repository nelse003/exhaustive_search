from __future__ import print_function

import logging
from collections import defaultdict

import iotbx
from giant.structure.restraints.occupancy import overlapping_occupancy_groups
from iotbx.pdb import hierarchy

logging = logging.getLogger(__name__)


def get_occupancy_group_grid_points(pdb, bound_states, ground_states, params):
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

    logging.info(
        "For all bound and ground states, "
        "select cartesian grid points for each altloc/residue \n"
        "involved in occupancy groups. A buffer of {} Angstrom \n"
        "is applied to minimal and maximal grid points,"
        "with a grid seperation of {}.\n".format(
            params.exhaustive.options.buffer, params.exhaustive.options.grid_spacing
        )
    )

    states = bound_states + ground_states

    pdb_in = iotbx.pdb.hierarchy.input(pdb)
    pdb_atoms = pdb_in.hierarchy.atoms()

    occupancy_group_cart_points = flex.vec3_double()
    for state in states:
        selection = state[0]
        selected_atoms = pdb_atoms.select(selection)
        sites_cart = selected_atoms.extract_xyz()
        grid_min = flex.double(
            [s - params.exhaustive.options.buffer for s in sites_cart.min()]
        )
        grid_max = flex.double(
            [s + params.exhaustive.options.buffer for s in sites_cart.max()]
        )
        grid_from_selection = grid.Grid(
            grid_spacing=params.exhaustive.options.grid_spacing,
            origin=tuple(grid_min),
            approx_max=tuple(grid_max),
        )

        logging.debug(grid_from_selection.summary())

        occupancy_group_cart_points = occupancy_group_cart_points.concatenate(
            grid_from_selection.cart_points()
        )

    logging.info(
        "Number of cartesian points to calculate "
        "|Fo-Fc| over: {}".format(len(occupancy_group_cart_points))
    )

    return occupancy_group_cart_points


def get_occupancy_groups(pdb, params):
    """
    Calculate occupancy groups given pdb file path.
    
    Wrapper of giant.structure.restraints.occupancy: overlapping_occupancy_groups(), 
    that generates hierarchy from pdb file path

    Parameters
    ----------

    :param pdb: 
    :param params: 

    Returns
    -------

    """

    logging.info("Gathering occupancy group information from PDB: %s", pdb)
    print("Gathering occupancy group information from PDB: %s", pdb)
    pdb_in = hierarchy.input(pdb)

    resnames = params.select.resnames.split(",")

    logging.info("Looking for ligands with resname {!s}".format(" or ".join(resnames)))

    occupancy_groups = overlapping_occupancy_groups(
        hierarchy=pdb_in.hierarchy,
        resnames=resnames,
        group_dist=params.select.group_dist,
        overlap_dist=params.select.overlap_dist,
        complete_groups=params.select.complete_groups,
        exclude_altlocs=params.select.exclude_altlocs.split(",")
        if params.select.exclude_altlocs
        else [],
        verbose=params.select.verbose,
    )

    return occupancy_groups


def get_selection_altloc_resseq_chain(sel_cache, altlocs, resseq, chain):

    """
    Get iotbx boolean selection of altlocs for a residue
    
    Given set of altlocs, residue number and chain,
    return a selection, using a supplied selection cache.

    Parameters
    ----------
    sel_cache: iotbx.pdb.atom_selection.cache
        selection cache on which to make the selection

    altlocs: list
        list of strings describing the altlocs involved i.e ['A','B']

    resseq: str
        residue sequnce number

    chain: str
        chain of supplied residue

    Returns
    --------
    altloc_selection, len(altlocs): tuple

        tuple containing the atoms in the ground state.
        altloc_selection is a iotbx.pdb selection object of type:

        scitbx_array_family_flex_ext.bool

    """

    selection_string = ""

    for altloc in altlocs:
        selection_string += "altloc {} or ".format(altloc)

    selection_string = selection_string.rstrip(" or ")

    if len(altlocs) > 1:
        selection_string = "({}) and chain {} and resseq {}".format(
            selection_string, chain, resseq
        )
    else:
        selection_string = "{} and chain {} and resseq {}".format(
            selection_string, chain, resseq
        )
    print(selection_string)

    altloc_selection = sel_cache.selection(selection_string)

    return (altloc_selection, len(altlocs))


# TODO Refactor into class
def get_bound_ground_states(pdb, params):
    """
    Get bound and ground states from PDB file.

    Uses the occupancy group information to get the ground and bound states

    Parameters
    -----------
    pdb: str
        path to pdb file

    params: libtbx.phil.scope_extract
            python object from phil file,
            edited with any additional parameters

    Returns
    -------
    ground_states: list
        list containing the atoms in the ground state.

        list is shaped:
        [[selection, number of altlocs],[selection, number of altlocs]...]

        where selection is a iotbx.pdb selection object of type:

        scitbx_array_family_flex_ext.bool

    bound_states: list
        list containing the atoms in the bound state.

        list is shaped:
        [[selection, number of altlocs],[selection, number of altlocs]...]

        where selection is a iotbx.pdb selection object of type:

        scitbx_array_family_flex_ext.bool

    Notes
    -----
    Currently the altlocs needs to be in multiple loops,
    such that altlocs are joint in selection,
    otherwise occupancy changes.
    I presume this might be to do with the
    distribution of points over which fo_fc is sampled?

    """
    logging.info("Process pdb file to get bound and ground states.")

    occupancy_groups = get_occupancy_groups(pdb, params)

    logging.debug(occupancy_groups)

    # Get bound altlocs if residue is in params.select.resname
    bound_altlocs = []
    for occupancy_group in occupancy_groups[0]:
        for residue_altloc in occupancy_group:

            if residue_altloc.get("resname") in params.select.resnames:
                bound_altlocs += residue_altloc.get("altloc")

    # Produce a default dictionary of residues
    # involved in ground and bound state,
    # and which altlocs belong to each state
    move_res = defaultdict(list)
    for occupancy_group in occupancy_groups[0]:
        for residue_altloc in occupancy_group:

            altloc = residue_altloc.get("altloc")
            chain = residue_altloc.get("chain")
            resseq = residue_altloc.get("resseq")

            if altloc in bound_altlocs:
                move_res[(chain, resseq, "Bound")].append(altloc)
            else:
                move_res[(chain, resseq, "Ground")].append(altloc)

    pdb_inp = iotbx.pdb.input(pdb)
    hier = pdb_inp.construct_hierarchy()
    sel_cache = hier.atom_selection_cache()

    bound_states = []
    ground_states = []
    for (chain, resseq, state), altlocs in move_res.iteritems():

        if state == "Bound":
            bound_states.append(
                get_selection_altloc_resseq_chain(sel_cache, altlocs, resseq, chain)
            )
        elif state == "Ground":
            ground_states.append(
                get_selection_altloc_resseq_chain(sel_cache, altlocs, resseq, chain)
            )

    try:
        ground_states
    except NameError:
        logging.info("There is no ground state. Try remodelling ground state")

    try:
        bound_states
    except NameError:
        logging.info("There is no bound state.")

    logging.info("len occ_groups {}".format(len(occupancy_groups)))

    logging.info("BOUND")
    logging.info(bound_states)

    logging.info("GROUND")
    logging.info(ground_states)

    return bound_states, ground_states
