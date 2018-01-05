from __future__ import print_function
import sys
import iotbx.pdb
from scitbx.array_family import flex

# Comparing bound & ground states
# TODO Add unit tests
# TODO Return atoms of interest rather than chains?

def select_chains_to_vary(bound_pdb_path, ground_pdb_path):

    """ Select difference between two pdb files"""

    # Get protein hierarchies
    bound_pdb_inp = iotbx.pdb.input(bound_pdb_path)
    ground_pdb_inp = iotbx.pdb.input(ground_pdb_path)
    bound_ph = bound_pdb_inp.construct_hierarchy()
    ground_ph = ground_pdb_inp.construct_hierarchy()

    # TODO add option to work with other atoms (not just hetatms)
    bound_hetatm_hier = select_all_hetatm(bound_ph)
    ground_hetatm_hier = select_all_hetatm(ground_ph)
    bound_no_water_hier = remove_waters(bound_hetatm_hier)
    ground_no_water_hier = remove_waters(ground_hetatm_hier)

    # TODO Add a try excpet block and excpetions to catch when there are not differnces in both bound and ground
    bound_not_in_ground_chains, ground_not_in_bound_chains = compare_hier(bound_no_water_hier,ground_no_water_hier)

    bound_not_in_ground_hier, sel_bound_not_in_ground_hier = chains_select_hier(bound_not_in_ground_chains, bound_ph)
    ground_not_in_bound_hier, sel_ground_not_in_bound_hier = chains_select_hier(bound_not_in_ground_chains, ground_ph)

    ## Are the differences between the bound and ground state within a rmsd cutoff
    if hierarchy_overlap(bound_not_in_ground_hier,  ground_not_in_bound_hier, rmsd_distance=5):
        ground_bound_chains = list(set(chain_ids_b + chain_ids_g))
    else:
        ground_bound_chains = []

    return ground_bound_chains, bound_not_in_ground_chains, ground_not_in_bound_chains


def select_all_hetatm(protein_hier, cache=None, copy=True):
    if not cache: cache = protein_hier.atom_selection_cache()
    sel =cache.selection('hetero')
    return protein_hier.select(sel, copy_atoms = copy)


def remove_waters(protein_hier, cache=None, copy=True):
    if not cache: cache=protein_hier.atom_selection_cache()
    sel = cache.selection('(not resname HOH)')
    return protein_hier.select(sel, copy_atoms=copy)


def compare_hier(hier_1, hier_2, cache_1 = None, cache_2 = None, rmsd_cutoff = 0.1):

    """Take two protein hierarchies and find unique chains and chains with atoms that differ between the hierarchies"""

    # Generate caches for each hierarchy
    if not cache_1: cache_1 = hier_1.atom_selection_cache()
    if not cache_2: cache_2 = hier_2.atom_selection_cache()

    # If hierarchy is the same return None
    if hier_1.is_similar_hierarchy(hier_2):
        print ("Hierarchies are the same")
        return [], []

    # Get chains from hierarchies
    chain_ids_1 = [chain.id for chain in hier_1.only_model().chains()]
    chain_ids_2 = [chain.id for chain in hier_2.only_model().chains()]

    # Find common, and unique chains in each hierarchy
    common_chain_ids = list(set(chain_ids_1).intersection(chain_ids_2))

    chain_id1_not_in_id2 = list(set(chain_ids_1) - set(chain_ids_2))
    chain_id2_not_in_id1 = list(set(chain_ids_2) - set(chain_ids_1))

    # Find chains with atoms that are sperated by more than rmsd cutoff
    chains_with_distant_atoms = []
    for chain_id in common_chain_ids:
        hier_1_sel = cache_1.selection("chain {}".format(chain_id))
        hier_2_sel = cache_2.selection("chain {}".format(chain_id))

        hier_1_xyz = hier_1.atoms().extract_xyz().select(hier_1_sel)
        hier_2_xyz = hier_2.atoms().extract_xyz().select(hier_2_sel)

        if len(hier_1_xyz) == len(hier_2_xyz):

            distance = flex.mean((hier_1_xyz-hier_2_xyz).dot())

            if distance > rmsd_cutoff**2:
                print("Chain {} is different by {} A RMSD".format(chain_id, distance))
                chains_with_distant_atoms.append(chain_id)
        else:
            print ("Chain {} is a different length in ground and bound states. "
                   "Something has been removed?. Assuming chain is different")
            chains_with_distant_atoms.append(chain_id)

    # Add chains that have atoms in different positions to chains that are not shared between hierarchies
    chains_hier_1_different_heir_2 = chains_with_distant_atoms + chain_id1_not_in_id2
    chains_hier_2_different_heir_1 = chains_with_distant_atoms + chain_id2_not_in_id1

    return chains_hier_1_different_heir_2, chains_hier_2_different_heir_1


def chains_select_hier(chain_list, hier, cache = None, copy = True):

    """Return hierarchy and selection including only selected chains"""

    if not cache: cache=hier.atom_selection_cache()

    chains_sel = []
    for chain in chain_list:
        chains_sel.append("chain {}".format(chain))

    chains_selection_str = ' and '.join(chains_sel)
    sel = cache.selection(chains_selection_str)

    return hier.select(sel, copy_atoms=copy), sel


def hierarchy_overlap(hier_1, hier_2, cache_1 = None, cache_2 = None, rmsd_distance=5):

    """ Returns True if any atom in two heirarchies are within given rmsd distance"""

    if not cache_1: cache_1 = hier_1.atom_selection_cache()
    if not cache_2: cache_2 = hier_2.atom_selection_cache()

    overlap = False
    for atom_hier_1 in hier_1.atoms():
        for atom_hier_2 in hier_2.atoms():

            distance = atom_hier_1.distance(atom_hier_2)

            if distance < rmsd_distance:
                overlap = True
                return overlap

    return overlap


def cmd_run(args, xtal_name):
    ground_pdb_path = "/dls/labxchem/data/2016/lb13385-64/processing/analysis_xce2_june_17/initial_model/DCP2B-x0146/Refine_0008/refine_8.split.ground-state.pdb"
    bound_pdb_path = "/dls/labxchem/data/2016/lb13385-64/processing/analysis_xce2_june_17/initial_model/DCP2B-x0146/Refine_0008/refine_8.split.bound-state.pdb"
    select_chains_to_vary(bound_pdb_path, ground_pdb_path)

if (__name__ == "__main__"):
    cmd_run(args=sys.argv[1:], xtal_name=None)