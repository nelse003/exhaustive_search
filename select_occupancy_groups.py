from __future__ import print_function
from iotbx.pdb import hierarchy
import iotbx
import libtbx.phil
import sys
from giant.structure.restraints.occupancy import overlapping_occupancy_groups
from giant.maths.geometry import pairwise_dists
from select_hierarchies import hierarchy_all_overlap
import itertools
from copy import deepcopy

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
# Shared Functions
##############################################################

# Process pdb file to provide occupancy groups
def get_occupancy_groups(pdb, params=master_phil.extract()):
    print("PDB",pdb)
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
    """ Powerset with a minimal length feature

    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    powerset([1,2,3],2) --> (1,2) (1,3) (2,3) (1,2,3)
    """

    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(min_len, len(s) + 1))

######################################################################
# Selection by residue  with conicident altlocs
######################################################################

def get_residue_altloc_dict(occupancy_groups):

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

    return residue_altloc_dict

def get_parameter_from_occupancy_groups(occupancy_groups, parameter_str):

    parameters = []

    for occupancy_group in occupancy_groups:
        for group in occupancy_group:
            for residue_altloc in group:
                if residue_altloc.get("model") == '':
                    parameters.append(residue_altloc.get(parameter_str))
                    #print("This parameter {} is not included in the occupancy group".format(parameter_str))
                else:
                    raise Warning("Multiple models are present in pdb file. "
                                  "This is not processable with occupancy group selection")

    return parameters

def within_rmsd_cutoff(atoms1, atoms2, params):
    for i in range(0, len(pairwise_dists(atoms1.extract_xyz(), atoms2.extract_xyz()))):
        if pairwise_dists(atoms1.extract_xyz(), atoms2.extract_xyz())[i][i] < params.settings.coincident_cutoff:
            continue
        else:
            return False
    return True

def get_altloc_groups_coincident(hier, residue_chain, altloc_group, params):

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

    for residue, altlocs in residue_altloc_dict.items():
        residue_powerset = []
        for altloc_group in powerset(altlocs, 2):
            residue_powerset.append(get_altloc_groups_coincident(hier, residue, altloc_group, params))
        yield residue_powerset


def get_state_selection(hier, coincident, occupancy_groups, params):

    lig_altloc_group = get_ligand_coincident_altloc_group(hier, coincident, params)
    sel_cache = hier.atom_selection_cache()

    coincident_loop = deepcopy(coincident)
    all_altlocs = set(get_parameter_from_occupancy_groups(occupancy_groups, "altloc"))

    bound_states = []
    ground_states = []
    for coincident_altloc_group in coincident_loop:
        if coincident_altloc_group[0] == lig_altloc_group:
            coincident.remove(coincident_altloc_group)
            print("Bound State", coincident_altloc_group)
            bound_states.append(get_bound_ground_selection(sel_cache,coincident_altloc_group))
            continue
        if list(all_altlocs.difference(lig_altloc_group)) == list(coincident_altloc_group[0]):
            coincident.remove(coincident_altloc_group)
            print("Ground State", coincident_altloc_group)
            ground_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))
            continue
    if coincident:
        print("\nThese altlocs do not follow the ground bound split")
        print(coincident)
        print("Using ligand altloc {} we presume these states to be bound\n".format(lig_altloc_group[0]))
        for coincident_altloc_group in coincident:
            if lig_altloc_group[0] in coincident_altloc_group[0]:
                print("Bound State", coincident_altloc_group)
                bound_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))
            else:
                print("Ground State", coincident_altloc_group)
                ground_states.append(get_bound_ground_selection(sel_cache, coincident_altloc_group))

    return bound_states, ground_states

def get_bound_ground_selection(sel_cache,coincident_altloc_group):

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

def process_refined_pdb_bound_ground_states(pdb, params=master_phil.extract()):

    occupancy_groups = get_occupancy_groups(pdb, params)
    pdb_inp = iotbx.pdb.input(pdb)
    hier = pdb_inp.construct_hierarchy()
    residue_altloc_dict = get_residue_altloc_dict(occupancy_groups)
    coincident = get_coincident_group(hier, residue_altloc_dict, params)
    bound_states, ground_states = get_state_selection(hier, coincident, occupancy_groups, params)

    return bound_states, ground_states

# def get_largest_concident_group(pdb, residue_altloc_dict, params):
#
#     coincident = get_coincident_group(pdb, residue_altloc_dict, params)
#
#     #print(residue_altloc_dict)
#
#     for residue_chain in residue_altloc_dict.keys():
#         residue = residue_chain[0]
#         chain = residue_chain[1]
#         altloc_group = residue_altloc_dict.get(residue_chain)
#         non_coincident_altlocs = list(altloc_group)
#         for conincident_altloc_group in coincident:
#             if conincident_altloc_group[1] == residue and conincident_altloc_group[2] == chain:
#                 for altloc in conincident_altloc_group[0]:
#                     non_coincident_altlocs.remove(altloc)
#         if non_coincident_altlocs:
#             raise ValueError("More than two groups of altlocs are present. This is currently not handled")

def get_ligand_coincident_altloc_group(hier, coincident, params):

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

######################################################################
# Selection using altlocs where all residues are shared, and they overlap
######################################################################

# Use occupancy groups to generate selection of atoms associated to bound & ground structure

# Return an atom selection on the protein hier to be used as selection
# Provide whether sufficent for co-varying

# Need to consider generation of grid (for summing FoFc) for seperated components.

###########################################################################
# Used in generation of cartesian point selection
###########################################################################
def generate_altloc_hiearchy(altloc_residue_dict, pdb):

    pdb_in = hierarchy.input(pdb)
    altloc_selection = generate_altloc_selection(altloc_residue_dict, pdb)
    altloc_hierarchy = pdb_in.hierarchy.select(altloc_selection)

    print_hier_atoms(altloc_hierarchy)

    return altloc_hierarchy

def generate_altloc_selection(altloc_residue_dict, pdb):

    pdb_in = hierarchy.input(pdb)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()

    selection_list = []
    for set in altloc_residue_dict.values():
        for resseq_id in set:
            selection_list.append("resid {} or ".format(resseq_id))

    selection_string = ''.join(selection_list)
    selection_string = selection_string.rstrip(" or ")
    selection_string = "altloc {} and ({})".format(str(altloc_residue_dict.keys()[0]), selection_string)

    altloc_selection = sel_cache.selection(selection_string)

    return altloc_selection

def from_altloc_generate_altloc_selection(pdb, altloc, occupancy_groups, params=master_phil.extract()):

    altloc_residue_dict = get_altloc_residue_dict(altloc, occupancy_groups)
    altloc_selection = generate_altloc_selection(altloc_residue_dict, pdb)

    return altloc_selection

def get_coincident_altloc_groups(altloc_groups_with_shared_residues, occupancy_groups, pdb):

    coincident_altloc_groups = []
    for altloc_group in altloc_groups_with_shared_residues:
        print("Altlocs with Shared residues {}".format(altloc_group))

        all_altloc_hier = []
        for altloc in altloc_group:
            altloc_residue_dict = get_altloc_residue_dict(altloc, occupancy_groups)
            altloc_hierarchy = generate_altloc_hiearchy(altloc_residue_dict, pdb)
            all_altloc_hier.append(altloc_hierarchy)

        # TODO Rewrite hierarchy overlap to handle any number of hierarchies

        if len(all_altloc_hier) > 2:
            raise Warning("Only two hierarchies can be compared currently. "
                          "Assuming that if a larger group, that they all overlap if the first two states overlap")

        if hierarchy_all_overlap(all_altloc_hier[0], all_altloc_hier[1], rmsd_distance=0.02):
            print("Altlocs {} are Coincident hierarchies".format(altloc_group))
            coincident_altloc_groups.append(altloc_group)

    return coincident_altloc_groups

def get_altloc_residue_dict(altlocs,occupancy_groups):

    altloc_residue_dict = dict()
    for altloc in altlocs:
        altloc_residue_set = set()
        for occupancy_group in occupancy_groups:
            for group in occupancy_group:
                for residue_altloc in group:
                    if altloc == residue_altloc.get("altloc"):
                        altloc_residue_set.add(residue_altloc.get("resseq"))
                        altloc_residue_dict.update({altloc: altloc_residue_set})

    return altloc_residue_dict


def find_altlocs_with_shared_residues(altlocs, occupancy_groups):

    altloc_residue_dict = get_altloc_residue_dict(altlocs,occupancy_groups)

    altloc_groups_with_shared_residues = list()
    for altloc_group in powerset(altloc_residue_dict.keys(), min_len=2):
        altloc_group_set = set()
        for altloc in altloc_group:
            altloc_group_set = altloc_group_set.union(altloc_residue_dict.get(altloc))

        different_residues_in_set = False
        for altloc in altloc_group:
            if altloc_group_set.difference(altloc_residue_dict.get(altloc)):
                different_residues_in_set = True
            else:
                pass
        if not different_residues_in_set:
            altloc_groups_with_shared_residues.append(altloc_group)

    return (altloc_groups_with_shared_residues)

# TODO Clean up how this sources the params
def get_coincident_altlocs(pdb, params=master_phil.extract()):

    occupancy_groups = get_occupancy_groups(pdb, params)
    altlocs = get_parameter_from_occupancy_groups(occupancy_groups, "altloc")
    altloc_groups_with_shared_residues = find_altlocs_with_shared_residues(altlocs, occupancy_groups)

    if len(altloc_groups_with_shared_residues) > 2:
        raise Warning("Current code can only handle two states")

    coincident_altloc_groups = get_coincident_altloc_groups(altloc_groups_with_shared_residues,
                                                            occupancy_groups, pdb)

    print ("Coincident_altloc_groups: {}".format(coincident_altloc_groups))

    return coincident_altloc_groups

# TODO Clean up how this sources the params
def get_altloc_hier(pdb, params=master_phil.extract()):
    coincident_altloc_groups = get_coincident_altlocs(pdb, params)
    occupancy_groups = get_occupancy_groups(pdb, params)

    for altloc_group in coincident_altloc_groups:
        altloc_residue_dict = get_altloc_residue_dict(altloc_group[0], occupancy_groups)
        altloc_hier = generate_altloc_hiearchy(altloc_residue_dict, pdb)

        yield altloc_hier

#######################################################################
# Helper functions
#######################################################################

def print_hier_atoms(hierarchy):

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
    pass



if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)