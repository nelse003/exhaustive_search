from __future__ import print_function

from copy import deepcopy

import iotbx
import itertools
from giant.maths.geometry import pairwise_dists
from giant.structure.restraints.occupancy import overlapping_occupancy_groups
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
    :type path
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


def powerset(iterable, min_len=0):
    """
    Powerset with a minimal length feature

    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    powerset([1,2,3],2) --> (1,2) (1,3) (2,3) (1,2,3)
    
    :param iterable: 
    :param min_len: 
    :return: 
    """
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(min_len, len(s) + 1))

######################################################################
# Selection by residue  with conicident altlocs
######################################################################

def get_residue_altloc_dict(occupancy_groups):
    """
    Return a dictionary of the format
    {(Residue, Chain):[altloc, altloc, altloc], ...} from occupancy groups
    
    :param occupancy_groups: 
    :return: 
    """

    residue_altloc_dict = {}

    residues = get_parameter_from_occupancy_groups(occupancy_groups, "resseq")
    chains = get_parameter_from_occupancy_groups(occupancy_groups, "chain")

    for residue_chain in zip(residues,chains):
        residue_altloc_set = set()
        for occupancy_group in occupancy_groups:
            for group in occupancy_group:
                for residue_altloc in group:
                    if residue_chain[0] == residue_altloc.get("resseq") \
                            and residue_chain[1] == residue_altloc.get("chain"):
                        residue_altloc_set.add(residue_altloc.get("altloc"))
                        residue_altloc_dict.update({residue_chain:
                                                        residue_altloc_set})

    logging.info("The residues defined in the occupancy group"
                " have atlocs:\n {}".format(residue_altloc_dict))

    return residue_altloc_dict

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

def get_altloc_groups_coincident(hier, residue_chain, altloc_group, params):
    """
    Determine whether all altlocs in an altloc group are coincident.
    
    :param hier: 
    :param residue_chain: 
    :param altloc_group: 
    :param params: 
    :return: 
    """

    sel_cache = hier.atom_selection_cache()
    pdb_atoms = hier.atoms()

    residue = residue_chain[0]
    chain = residue_chain[1]

    for altloc, compare_altloc in itertools.combinations(altloc_group, 2):

        selection_string = "altloc {} and resseq {} and chain {}".format(altloc, residue, chain)
        altloc_selection = sel_cache.selection(selection_string)
        altloc_atoms = pdb_atoms.select(altloc_selection)

        compare_selection_string = "altloc {} and resseq {} and chain {}".format(compare_altloc, residue, chain)
        compare_altloc_selection = sel_cache.selection(compare_selection_string)
        compare_altloc_atoms = pdb_atoms.select(compare_altloc_selection)

        if within_rmsd_cutoff(altloc_atoms,compare_altloc_atoms, params):
            continue
        else:
            return (altloc_group, residue, chain, False)
    return (altloc_group, residue, chain, True)

def get_coincident_residue_powersets(hier, residue_altloc_dict, params):
    """
    For altlocs of each residue determine whether a powerset of altlocs are coincident.
    
    Given a residue with altlocs [A,B,C,D]. i.e {(121, A): [A, B, C, D]}, yield a list of altloc groups 
    that are coincident:
    
    [ ([A,B,C,D], 121, A, False), 
      ([A,B,C], 121, A, False),
      ([A,C,D], 121, A, False),
      ([B,C,D], 121, A, False),
      ([B,A,D], 121, A, False),
      ([A,B], 121, A, True),
      ([C,D], 121, A, True),... ]
      
    :param hier: 
    :param residue_altloc_dict: 
    :type residue_altloc_dict: dict
    :param params: 
    :return: 
    """

    for residue, altlocs in residue_altloc_dict.items():
        residue_powerset = []
        for altloc_group in powerset(altlocs, 2):
            residue_powerset.append(get_altloc_groups_coincident(hier, residue, altloc_group, params))
        yield residue_powerset


def get_largest_coincident(coincident,occupancy_groups):

    """ Take a powerset of coincident altlocs with residue and chain. Return the largest set(s) of altlocs.
     
     Will return single set if contains all altlocs for that resiude, or multiple if there are multiple groups.
    
    :param coincident: 
    :param occupancy_groups: 
    :return: 
    """

    residues = get_parameter_from_occupancy_groups(occupancy_groups, "resseq")
    chains = get_parameter_from_occupancy_groups(occupancy_groups, "chain")

    residue_chains = set()
    for residue_chain in zip(residues, chains):
        residue_chains.add(residue_chain)

    largest_coincident = []
    for residue_chain in residue_chains:

        list_all_coincident_altlocs = []
        all_altocs_res_chain = set()
        for coincident_altlocs in coincident:
            if coincident_altlocs[1] == residue_chain[0] and coincident_altlocs[2] == residue_chain[1]:
                list_all_coincident_altlocs.append(set(coincident_altlocs[0]))
                for altloc in coincident_altlocs[0]:
                    all_altocs_res_chain.add(altloc)

        list_all_coincident_altlocs.sort(key=len, reverse=True)
        for altloc_set in list_all_coincident_altlocs:
            difference_set = altloc_set.symmetric_difference(all_altocs_res_chain)
            if difference_set == set():
                largest_coincident.append((tuple(altloc_set),residue_chain[0], residue_chain[1]))
                break
            elif difference_set in list_all_coincident_altlocs:
                largest_coincident.append((tuple(altloc_set),residue_chain[0], residue_chain[1]))
                largest_coincident.append((tuple(difference_set), residue_chain[0], residue_chain[1]))
                break
            else:
                logging.info("The set of altlocs do not correspond to a two state system. Trying a shorter set of altlocs")

    return largest_coincident


def get_state_selection(hier, coincident, occupancy_groups, params):
    """
    Determine a split between ground and bound states, given list of largest coincident altloc groups
    
    :param hier: 
    :param coincident: 
    :param occupancy_groups: 
    :param params: 
    :return: 
    """

    lig_altloc_group = get_ligand_coincident_altloc_group(hier, coincident, params)
    sel_cache = hier.atom_selection_cache()

    coincident_loop = deepcopy(coincident)
    all_altlocs = set(get_parameter_from_occupancy_groups(occupancy_groups, "altloc"))

    logging.info("Given coincident altloc groups return a split between ground and bound state \n" 
                "For altloc groups that follow the ligand i.e LIG altloc C and LIG altloc D \n" 
                " are coincident, then if RES 121 [(A,B) ,(C,D)], set (C,D) as bound state and \n" 
                "(A,B) as ground state \n")

    bound_states = []
    ground_states = []

    logging.debug("Lig altloc group: {}".format(lig_altloc_group))

    if coincident:
        for coincident_altloc_group in coincident_loop:

            logging.debug(coincident_altloc_group)

            if coincident_altloc_group[0] == lig_altloc_group:
                coincident.remove(coincident_altloc_group)
                logging.info("Bound State: {}".format(coincident_altloc_group))
                bound_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))
                continue
            if list(all_altlocs.difference(lig_altloc_group)) == list(coincident_altloc_group[0]):
                coincident.remove(coincident_altloc_group)
                logging.info("Ground State: {}".format(coincident_altloc_group))
                ground_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))
                continue
        if coincident:
            logging.info("These altlocs do not follow the ground bound split: \n {}\n".format(coincident) + \
                        "Using ligand altloc {} we presume these states to be bound\n".format(lig_altloc_group[0]))

            for coincident_altloc_group in coincident:
                if lig_altloc_group[0] in coincident_altloc_group[0]:
                    logging.info("Bound State: {}".format(coincident_altloc_group))
                    bound_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))
                else:
                    logging.info("Ground State: {}".format(coincident_altloc_group))
                    ground_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))

    if not ground_states:
        logging.warning("No Ground State detected")
    if not bound_states:
        logging.warning("No Bound State detected")

    logging.info("____________________________________________________________________________________________________")

    return bound_states, ground_states


def get_bound_ground_selection(sel_cache, coincident_altloc_group):
    """
    Get iotbx boolean selection of altlocs for a residue
    
    Given a coincident altloc group in format (('C', 'D'), ' 121', 'A') return a selection, 
    using a supplied selection cache.
    Also returns number of altlocs in that selection. In the format: [altloc_selection, num_altlocs]
    
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

def get_coincident_group(hier, residue_altloc_dict, params):
    """
    Get a list of coincident altloc groups.
    
    In format [(altloc_group,residue,chain),(altloc_group,residue,chain)...] 
    
    :param hier: 
    :param residue_altloc_dict: 
    :param params: 
    :return: 
    """

    logging.info("Get altlocs groups that are coincident, in format" 
                "[(altloc_group,residue,chain),(altloc_group,residue,chain)...] ")

    coincident = []
    for residue_powerset in get_coincident_residue_powersets(hier, residue_altloc_dict, params):
        for altloc_groups_coincident in residue_powerset:
            altloc_group = altloc_groups_coincident[0]
            residue = altloc_groups_coincident[1]
            chain = altloc_groups_coincident[2]
            coincident_flag = altloc_groups_coincident[3]
            if coincident_flag:
                coincident.append((altloc_group, residue, chain))
    return coincident


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

    print("HERE")

    logging.info("Process pdb file to get bound and ground states.")

    occupancy_groups = get_occupancy_groups(pdb, params)

    logging.debug(occupancy_groups)

    pdb_inp = iotbx.pdb.input(pdb)
    hier = pdb_inp.construct_hierarchy()
    sel_cache = hier.atom_selection_cache()

    if len(occupancy_groups) == 1 and len(occupancy_groups[0]) >= 2:

        print("THERE")

        # There are no coincident residues therefore use the occupancy id instead.
        # This may need doing earlier, i.e the selection by altlocs is possibly really stupid.
        # Only the case for if a single complete group (i.e first part of the if statement)

        # (altloc_group, residue, chain)
        # (('C', 'D'), ' 121', 'A')
        #[altloc_selection, num_altlocs]
        #{'chain': 'A', 'altloc': 'A', 'resseq': '  67', 'icode': ' ', 'resname': 'ARG', 'model': ''}

        bound_states = []
        ground_states = []
        move_res = dict()
        for occupancy_group in occupancy_groups[0]:

            bound_state_flag = False
            state = []
            for residue_altloc in occupancy_group:
                if residue_altloc.get('resname') in params.select.resnames:
                    bound_state_flag = True
                    state_string = "Bound"
                else:
                    state_string = "Ground"

            for residue_altloc in occupancy_group:

                altloc = residue_altloc.get('altloc')
                chain = residue_altloc.get('chain')
                resseq = residue_altloc.get('resseq')

                if move_res.has_key((chain,resseq)):
                    move_res[(chain,resseq)].append(altloc)
                else:
                    move_res[(chain,resseq)] = [altloc]

        print("XX {}".format(move_res))

        for residue_chain, altlocs in d.iteritems():

            print(residue_chain, altlocs)
            continue

            altloc = residue_altloc.get('altloc')
            chain = residue_altloc.get('chain')
            resseq = residue_altloc.get('resseq')
            logging.info("{} State: {}".format(state_string, (((altloc,), resseq, chain))))
            state.append(get_bound_ground_selection(sel_cache, (((altloc,), resseq, chain))))
            logging.debug("APPEND STATE")
            logging.debug(state)

        if bound_state_flag:
            bound_states += state
        else:
            ground_states += state

        exit()

        # try:
        #     ground_states
        # except NameError:
        #     logging.info("There is no ground state. Try remodelling ground state")
        # try:
        #     bound_states
        # except NameError:
        #     logging.info("There is no bound state.")

        logging.info("len occ_groups {}".format(len(occupancy_groups)))

        logging.info("BOUND")
        logging.info(bound_states)

        logging.info("GROUND")
        logging.info(ground_states)

        return bound_states, ground_states

    #TODO Add more appropritate elif to catch other more complex cases (> two states/ other ways to get to two states) #63
    else:

        print("ACTUALLY this is probably a untested piece of code, only come here if desperate")
        exit()

        residue_altloc_dict = get_residue_altloc_dict(occupancy_groups)
        coincident = get_coincident_group(hier, residue_altloc_dict, params)
        logging.debug("Coincident: {}".format(coincident))
        largest_coincident = get_largest_coincident(coincident,occupancy_groups)
        bound_states, ground_states = get_state_selection(hier, largest_coincident, occupancy_groups, params)
        return bound_states, ground_states


#######################################################################
# Helper functions
#######################################################################

def print_hier_atoms(hierarchy):
    """
    Basic printing of an iotbx.pdb.hierarchy
    
    :param hierarchy: 
    :return: 
    """

    for model in hierarchy.models():
        print("Model: {}".format(model.id))

        for chain in model.chains():
            print("Chain: {}".format(chain.id))

            for residue_group in chain.residue_groups():
                print("Residue: {}".format(residue_group.resseq))

                for atom_group in residue_group.atom_groups():
                    print("Altloc: {}".format(atom_group.altloc))

                    for atom in atom_group.atoms():
                        print("Atom Name: {}".format(atom.name))

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
