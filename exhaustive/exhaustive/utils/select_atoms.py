from __future__ import print_function

from copy import deepcopy

import itertools

from giant.maths.geometry import pairwise_dists
from giant.structure.restraints.occupancy import overlapping_occupancy_groups
import giant.grid as grid

import iotbx
from scitbx.array_family import flex
from iotbx.pdb import hierarchy
##############################################################
import logging
logging = logging.getLogger(__name__)
##############################################################
# Shared Functions
##############################################################

# Process pdb file to provide occupancy groups
def get_occupancy_groups(pdb, params):
    """
    Calculate occupancy groups given pdb file path.
    
    Wrapper of giant.structure.restraints.occupancy: overlapping_occupancy_groups(), 
    that generates hierarchy from pdb file path
    
    :param pdb: 
    :param params: 
    :return: 
    """

    logging.info("Gathering occupancy group information from PDB: %s", pdb)
    print("Gathering occupancy group information from PDB: %s", pdb)
    pdb_in = hierarchy.input(pdb)

    resnames = params.select.resnames.split(',')

    logging.info('Looking for ligands with resname {!s}'.format(' or '.join(resnames)))

    occupancy_groups = overlapping_occupancy_groups(hierarchy=pdb_in.hierarchy,
                                                    resnames=resnames,
                                                    group_dist=params.select.group_dist,
                                                    overlap_dist=params.select.overlap_dist,
                                                    complete_groups=params.select.complete_groups,
                                                    exclude_altlocs=params.select.exclude_altlocs.split(
                                                        ',') if params.select.exclude_altlocs else [],
                                                    verbose=params.select.verbose)

    return occupancy_groups


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


def get_bound_ground_selection(sel_cache, coincident_altloc_group):
    """
    Get iotbx boolean selection of altlocs for a residue
    
    Given a coincident altloc group in format (('C', 'D'), ' 121', 'A') return a selection, 
    using a supplied selection cache.
    Also returns number of altlocs in that selection. 
    In the format: [altloc_selection, num_altlocs]
    
    :param sel_cache: 
    :param coincident_altloc_group: 
    :return: 
    """

    print(coincident_altloc_group)
    print(coincident_altloc_group[0])

    num_altlocs =  len(coincident_altloc_group[0])
    selection_string = ""
    for altloc in coincident_altloc_group[0]:
        selection_string += "altloc {} or ".format(altloc)
    selection_string = selection_string.rstrip(" or ")

    if num_altlocs > 1:
        selection_string = "({}) and chain {} and resseq {}".format(selection_string, coincident_altloc_group[2],
                                                            coincident_altloc_group[1])
    else:
        selection_string = "{} and chain {} and resseq {}".format(selection_string, coincident_altloc_group[2], coincident_altloc_group[1])

    print(selection_string)
    altloc_selection = sel_cache.selection(selection_string)

    return [altloc_selection, num_altlocs]


def get_ligand_coincident_altloc_group(hier, coincident, params):
    """
    Determine whether altloc group is bound or ground. 
    
    Returns altloc group i.e (C,D) that has the same altloc as the first altloc of the ligand in the specified 
    ligand resnames (DRG,FRG,LIG,UNK,UNL). Used in cases where:
    
    LIG (C,D)
    Ground (A,B)
    
    Residue (A,C) (B,D)
    
    to pick the bound state of the residue@
    
    Bound (A,C)
    Ground (B,D)
    
    :param hier: 
    :param coincident: 
    :param params: 
    :return: 
    """

    resnames = params.select.resnames.split(',')
    sel_cache = hier.atom_selection_cache()
    pdb_atoms = hier.atoms()

    logging.debug("Resnames: {}".format(resnames))

    lig_altlocs = set()
    lig_chain = set()
    for resname in resnames:
        selection_string = "resname {}".format(resname)
        altloc_selection = sel_cache.selection(selection_string)
        altloc_atoms = pdb_atoms.select(altloc_selection)
        for atom in altloc_atoms:
            lig_altlocs.add(atom.parent().altloc)
            lig_chain.add(atom.parent().parent().parent().id)

    logging.debug("Lig altlocs: {}".format(lig_altlocs))
    logging.debug("Lig chain: {}".format(lig_chain))

    logging.debug("coincident: {}".format(coincident))

    for coincident_altloc_group in coincident:

        logging.debug("Coinicident altloc group: {}".format(coincident_altloc_group))

        print(coincident_altloc_group[2])
        print(list(lig_chain)[0])
        print(list(lig_altlocs)[0])
        print(coincident_altloc_group[0])

        if coincident_altloc_group[2] == list(lig_chain)[0] and list(lig_altlocs)[0] in coincident_altloc_group[0]:
            return coincident_altloc_group[0]
        else:
            continue
    # If no matching coinicdent cases have been made return just the ligand altlocs
    return tuple(lig_altlocs)


def process_refined_pdb_bound_ground_states(pdb, params):
    """
    Main Function that returns bound and ground states from occupancy group of the supplied PDB file.
    
    Returns in the form [[selection, number of altlocs],[selection, number of altlocs]...] 
    for each of the bound and ground states. 
    
    :param pdb: 
    :param params: 
    :return: 
    """
    logging.info("Process pdb file to get bound and ground states.")

    occupancy_groups = get_occupancy_groups(pdb, params)

    logging.debug(occupancy_groups)

    pdb_inp = iotbx.pdb.input(pdb)
    hier = pdb_inp.construct_hierarchy()
    sel_cache = hier.atom_selection_cache()

    if len(occupancy_groups) == 1 and len(occupancy_groups[0]) >= 2:

        # There are no coincident residues therefore use the occupancy id instead.
        # This may need doing earlier, i.e the selection by altlocs is possibly really stupid.
        # Only the case for if a single complete group (i.e first part of the if statement)

        # (altloc_group, residue, chain)
        # (('C', 'D'), ' 121', 'A')
        # [altloc_selection, num_altlocs]
        # {'chain': 'A', 'altloc': 'A', 'resseq': '  67', 'icode': ' ', 'resname': 'ARG', 'model': ''}

        bound_states = []
        ground_states = []
        move_res = dict()
        bound_altlocs = []
        for occupancy_group in occupancy_groups[0]:
            for residue_altloc in occupancy_group:

                if residue_altloc.get('resname') in params.select.resnames:
                    bound_altlocs += residue_altloc.get('altloc')

        for occupancy_group in occupancy_groups[0]:
            for residue_altloc in occupancy_group:

                altloc = residue_altloc.get('altloc')
                chain = residue_altloc.get('chain')
                resseq = residue_altloc.get('resseq')

                if altloc in bound_altlocs:
                    state_string = "Bound"
                else:
                    state_string = "Ground"

                if move_res.has_key((chain, resseq, state_string)):
                    move_res[(chain, resseq, state_string)].append(altloc)
                else:
                    move_res[(chain, resseq, state_string)] = [altloc]

        print("move_res: {}".format(move_res))

        for residue_chain, altlocs in move_res.iteritems():

            resseq = residue_chain[1]
            chain = residue_chain[0]
            state = residue_chain[2]

            logging.info("{} State: {}".format(state, ((tuple(altlocs), resseq, chain))))

            if state == "Bound":
                bound_states.append(get_bound_ground_selection(sel_cache, ((tuple(altlocs), resseq, chain))))
            elif state == "Ground":
                ground_states.append(get_bound_ground_selection(sel_cache, ((tuple(altlocs), resseq, chain))))
            else:
                raise ValueError("{} states is undefined".format(state))
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
