from __future__ import print_function

import sys
from copy import deepcopy

import iotbx
import itertools
import libtbx.phil
from giant.maths.geometry import pairwise_dists
from giant.structure.restraints.occupancy import overlapping_occupancy_groups
from iotbx.pdb import hierarchy

blank_arg_prepend = {'.pdb': 'pdb=', '.mtz': 'mtz=', '.csv': 'csv='}
##############################################################
master_phil = libtbx.phil.parse("""
input{
    pdb = None
        .type = path
    mtz = None
        .type = path
    csv = None
        .type = path
}
occupancy{
    resnames = DRG,FRG,LIG,UNK,UNL
        .help = 'Residues to generate constraint groups around for occupancy refinement (comma separated list of residue identifiers, i.e. resname=LIG or resname=LIG,UNL)'
        .type = str

    group_dist = 5
        .type = float
        .help = 'Distance to use when clustering atoms that should have the SAME occupancy'

    overlap_dist = 2
        .type = float
        .help = 'Distance to use when clustering atoms that should have occupancies that SUM TO LESS THAN ONE'

    exclude_altlocs = None
        .help = 'Exclude certain altlocs from occupancy groups (e.g. A or A,B)'
        .type = str

    complete_groups = True
        .help = 'Generate a set of fully constrained groups (that sum to unitary occupancy) when True. Generate a set of weaker constraints for overlapping atoms when False.'
        .type = bool
}        
output{
    out_dir = "/hdlocal/home/enelson/Dropbox/DPhil/exhaustive_search/output/"
        .type = str
}        
settings{
    verbose = True
        .type = bool
    coincident_cutoff = 0.05
        .help = ' RMSD Cutoff in Angstrom, for two structures considered coincident'
        .type = float
}

""", process_includes=True)
##############################################################
import logging

logger = logging.getLogger(__name__)
##############################################################
# Shared Functions
##############################################################

# Process pdb file to provide occupancy groups
def get_occupancy_groups(pdb, params=master_phil.extract()):
    """
    Calculate occupancy groups given pdb file path.
    
    Wrapper of giant.structure.restraints.occupancy: overlapping_occupancy_groups(), 
    that generates hierarchy from pdb file path
    
    :param pdb: 
    :type path
    :param params: 
    :return: 
    """

    logger.info("Gathering occupancy group information from PDB", pdb)
    pdb_in = hierarchy.input(pdb)

    resnames = params.occupancy.resnames.split(',')
    # TODO Change to log, not print statements
    if params.settings.verbose:
        print('Looking for ligands with resname {!s}'.format(' or '.join(resnames)))

    occupancy_groups = overlapping_occupancy_groups(hierarchy=pdb_in.hierarchy,
                                                    resnames=resnames,
                                                    group_dist=params.occupancy.group_dist,
                                                    overlap_dist=params.occupancy.overlap_dist,
                                                    complete_groups=params.occupancy.complete_groups,
                                                    exclude_altlocs=params.occupancy.exclude_altlocs.split(
                                                        ',') if params.occupancy.exclude_altlocs else [],
                                                    verbose=params.settings.verbose)

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
    Return a dictionary of the format {(Residue, Chain):[altloc, altloc, altloc], ...} from occupancy groups
    
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
                    if residue_chain[0] == residue_altloc.get("resseq") and residue_chain[1] == residue_altloc.get("chain"):
                        residue_altloc_set.add(residue_altloc.get("altloc"))
                        residue_altloc_dict.update({residue_chain: residue_altloc_set})

    logger.info("The residues defined in the occupancy group have atlocs:\n {}".format(residue_altloc_dict))

    return residue_altloc_dict

def get_parameter_from_occupancy_groups(occupancy_groups, parameter_str):
    """
    Extracts a parameter from occupancy groups, given a str matching that parameter (altloc, chain, resseq...)
    
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
                                  "This is not processable with occupancy group selection")
    if parameters == []:
        raise Warning("Parameter may not be recognised,as output list is empty")

    return parameters

def within_rmsd_cutoff(atoms1, atoms2, params):
    """
    Given two groups of atoms determine within a given cutoff (supplied via params)
    
    :param atoms1: 
    :param atoms2: 
    :param params: 
    :return: 
    """

    for i in range(0, len(pairwise_dists(atoms1.extract_xyz(), atoms2.extract_xyz()))):
        if pairwise_dists(atoms1.extract_xyz(), atoms2.extract_xyz())[i][i] < params.settings.coincident_cutoff:
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

        list_all_coincident_altlocs.sort(key = len, reverse =True)
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
                logger.info("The set of altlocs do not correspond to a two state system. Trying a shorter set of altlocs")

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

    logger.info("Given coincident altloc groups return a split between ground and bound state \n" \
                "For altloc groups that follow the ligand i.e LIG altloc C and LIG altloc D \n" \
                " are coincident, then if RES 121 [(A,B) ,(C,D)], set (C,D) as bound state and \n" \
                "(A,B) as ground state \n")

    bound_states = []
    ground_states = []

    print("lig_alt:",lig_altloc_group)

    for coincident_altloc_group in coincident_loop:

        print(coincident_altloc_group)

        if coincident_altloc_group[0] == lig_altloc_group:
            coincident.remove(coincident_altloc_group)
            logger.info("Bound State: {}".format(coincident_altloc_group))
            bound_states.append(get_bound_ground_selection(sel_cache,coincident_altloc_group))
            continue
        if list(all_altlocs.difference(lig_altloc_group)) == list(coincident_altloc_group[0]):
            coincident.remove(coincident_altloc_group)
            logger.info("Ground State: {}".format(coincident_altloc_group))
            ground_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))
            continue
    if coincident:
        logger.info("These altlocs do not follow the ground bound split: \n {}\n".format(coincident) + \
                    "Using ligand altloc {} we presume these states to be bound\n".format(lig_altloc_group[0]))

        for coincident_altloc_group in coincident:
            if lig_altloc_group[0] in coincident_altloc_group[0]:
                logger.info("Bound State: {}".format(coincident_altloc_group))
                bound_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))
            else:
                logger.info("Ground State: {}".format(coincident_altloc_group))
                ground_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))

    return bound_states, ground_states

def get_bound_ground_selection(sel_cache,coincident_altloc_group):
    """
    Get iotbx boolean selection of altlocs for a residue
    
    Given a coincident altloc group in format (('C', 'D'), ' 121', 'A') return a selection, 
    using a supplied selection cache.
    Also returns number of altlocs in that selection.
    
    :param sel_cache: 
    :param coincident_altloc_group: 
    :return: 
    """

    num_altlocs =  len(coincident_altloc_group[0])
    selection_string = ""
    for altloc in coincident_altloc_group[0]:
        selection_string += "altloc {} or ".format(altloc)
    selection_string = selection_string.rstrip(" or ")
    selection_string = "({}) and chain {} and resseq {}".format(selection_string, coincident_altloc_group[2],
                                                                coincident_altloc_group[1])
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

    logger.info("Get altlocs groups that are coincident, in format" \
                "[(altloc_group,residue,chain),(altloc_group,residue,chain)...] ")

    coincident = []
    for residue_powerset in get_coincident_residue_powersets(hier, residue_altloc_dict, params):
        for altloc_groups_coincident in residue_powerset:
            altloc_group = altloc_groups_coincident[0]
            residue = altloc_groups_coincident[1]
            chain = altloc_groups_coincident[2]
            coincident_flag = altloc_groups_coincident[3]
            if coincident_flag:
                coincident.append((altloc_group,residue,chain))
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

    resnames = params.occupancy.resnames.split(',')
    sel_cache = hier.atom_selection_cache()
    pdb_atoms = hier.atoms()

    lig_altlocs = set()
    lig_chain = set()
    for resname in resnames:
        selection_string = "resname {}".format(resname)
        altloc_selection = sel_cache.selection(selection_string)
        altloc_atoms = pdb_atoms.select(altloc_selection)
        for atom in altloc_atoms:
            lig_altlocs.add(atom.parent().altloc)
            lig_chain.add(atom.parent().parent().parent().id)

    for coincident_altloc_group in coincident:

        if coincident_altloc_group[2] == list(lig_chain)[0] and list(lig_altlocs)[0] in coincident_altloc_group[0]:
            return coincident_altloc_group[0]
        else:
            continue

def process_refined_pdb_bound_ground_states(pdb, params=master_phil.extract()):
    """
    Main Function that returns bound and ground states from occupancy group of the supplied PDB file.
    
    Returns in the form [[selection, number of altlocs],[selection, number of altlocs]...] 
    for each of the bound and ground states. 
    
    :param pdb: 
    :param params: 
    :return: 
    """

    logger.info("Process pdb file to get bound and ground states.")

    occupancy_groups = get_occupancy_groups(pdb, params)
    pdb_inp = iotbx.pdb.input(pdb)
    hier = pdb_inp.construct_hierarchy()
    residue_altloc_dict = get_residue_altloc_dict(occupancy_groups)
    coincident = get_coincident_group(hier, residue_altloc_dict, params)
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

    for chain in hierarchy.only_model().chains():
        print("Chain: {}".format(chain.id))

        for residue_group in chain.residue_groups():
            print("Residue: {}".format(residue_group.resseq))

            for atom_group in residue_group.atom_groups():
                print("Altloc: {}".format(atom_group.altloc))

                for atom in atom_group.atoms():
                    print("Atom Name: {}".format(atom.name))

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

def run(params):
    occupancy_groups = get_occupancy_groups(params.input.pdb)
    get_parameter_from_occupancy_groups(occupancy_groups, "CAT")


if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)