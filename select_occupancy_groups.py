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