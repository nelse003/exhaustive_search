from __future__ import print_function
import sys
import iotbx.pdb

### Comparing bound & ground states ######

def select_atoms_to_vary(bound_pdb_path, ground_pdb_path):
    # Get protein hierarchies
    bound_pdb_inp = iotbx.pdb.input(bound_pdb_path)
    ground_pdb_inp = iotbx.pdb.input(ground_pdb_path)
    bound_ph = bound_pdb_inp.construct_hierarchy()
    ground_ph = ground_pdb_inp.construct_hierarchy()

    bound_hetatm_hier = select_all_hetatm(bound_ph)
    ground_hetatm_hier = select_all_hetatm(ground_ph)
    bound_no_water_hier = remove_waters(bound_hetatm_hier)
    ground_no_water_hier = remove_waters(ground_hetatm_hier)

    print(ground_no_water_hier.overall_counts().as_str())
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(bound_no_water_hier.overall_counts().as_str())

    bound_not_in_ground_hier, ground_not_in_bound_hier = compare_hier(bound_no_water_hier,
                                                                       ground_no_water_hier)
    # if ground_not_in_bound_hier is not None:
    #     raise ValueError("Ground State has atoms not in bound state")
    #
    # if bound_not_in_ground_hier is None:
    #     raise ValueError("Bound state has no differing atoms")
    #
    # near_bound_differences_hier = select_nearby_hier(bound_not_in_ground_hier, distance=5, hetatm_only=True)
    #
    # if near_bound_differences_hier is None:
    #     is_bound_only = True
    #     asc = bound_no_water_hier.atom_selection_cache()
    # else:
    #     is_bound_only = False
    #     join_two_hier(near_bound_differences_hier, bound_no_water_hier)
    #     asc = near_bound_differences_hier.atom_selection_cache()
    #
    # sel_atoms = asc.selection("all")
    #
    # return sel_atoms, is_bound_only


def select_all_hetatm(protein_hier, cache=None, copy=True):
    if not cache: cache = protein_hier.atom_selection_cache()
    sel =cache.selection('hetero')
    return protein_hier.select(sel, copy_atoms = copy)


def remove_waters(protein_hier, cache=None, copy=True):
    if not cache: cache=protein_hier.atom_selection_cache()
    sel = cache.selection('(not resname HOH)')
    return protein_hier.select(sel, copy_atoms=copy)


def compare_hier(hier_1, hier_2, cache_1 = None, cache_2 = None):

    if not cache_1: cache_1 = hier_1.atom_selection_cache()
    if not cache_2: cache_2 = hier_2.atom_selection_cache()

    if hier_1.is_similar_hierarchy(hier_2):
        return None, None

    for chain in hier_1.only_model().chains():
        print(chain.id)

    return heir_1_atoms_not_in_hier_2, heir_2_atoms_not_in_hier_1


def join_two_hier(protein_hier_1, protein_hier_2):
    return joined_hier


def select_nearby_hier(protein_hier, distance=5, hetatm_only=False):
    return nearby_hier


def compare_atom_positions(selection_ground, selection_bound, rmsd_cutoff):

    return atoms_in_different_positions


def all_atoms_in_res(protein_hier):

    return all_atom_hier

def cmd_run(args, xtal_name):
    ground_pdb_path = "/dls/labxchem/data/2016/lb13385-64/processing/analysis_xce2_june_17/initial_model/DCP2B-x0146/Refine_0008/refine_8.split.ground-state.pdb"
    bound_pdb_path = "/dls/labxchem/data/2016/lb13385-64/processing/analysis_xce2_june_17/initial_model/DCP2B-x0146/Refine_0008/refine_8.split.bound-state.pdb"
    select_atoms_to_vary(bound_pdb_path, ground_pdb_path)

if (__name__ == "__main__"):
    cmd_run(args=sys.argv[1:], xtal_name=None)