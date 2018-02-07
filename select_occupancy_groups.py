from __future__ import print_function
from iotbx.pdb import hierarchy
import libtbx.phil
import sys
from giant.structure.restraints.occupancy import overlapping_occupancy_groups
from giant.structure.altlocs import find_duplicate_conformers
from select_hierarchies import hierarchy_all_overlap
import itertools

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
}

""", process_includes=True)


##############################################################


# Process pdb file to provide occupancy groups
def get_occupancy_groups(pdb, params):
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


# Use occupancy groups to generate selection of atoms associated to bound & ground structure

# Return an atom selection on the protein hier to be used as selection
# Provide whether sufficent for co-varying

# Need to consider generation of grid (for summing FoFc) for seperated components.

def get_altlocs_from_occupancy_groups(occupancy_groups):

    altlocs = set()

    for occupancy_group in occupancy_groups:
        for group in occupancy_group:
            for residue_altloc in group:
                if residue_altloc.get("model") == '':
                    altlocs.add(residue_altloc.get("altloc"))
                else:
                    raise Warning("Multiple models are present in pdb file. "
                                  "This is not processable with occupancy group selection")

    return list(altlocs)

def powerset(iterable,min_len=0):
    """ Powerset with a minimal length feature
    
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    powerset([1,2,3],2) --> (1,2) (1,3) (2,3) (1,2,3)
    """

    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(min_len,len(s)+1))

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

def generate_altloc_hiearchy(altloc_residue_dict, pdb):

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
    altloc_hierarchy = pdb_in.hierarchy.select(altloc_selection)

    return altloc_hierarchy

def get_coincident_altloc_groups(altloc_groups_with_shared_residues, occupancy_groups, pdb):

    coincident_altloc_groups = []
    for altloc_group in altloc_groups_with_shared_residues:
        all_altloc_hier = []
        for altloc in altloc_group:
            altloc_residue_dict = get_altloc_residue_dict(altloc, occupancy_groups)
            altloc_hierarchy = generate_altloc_hiearchy(altloc_residue_dict, pdb)
            all_altloc_hier.append(altloc_hierarchy)

        # TODO Rewrite hierarchy overlap to handle any number of hierarchies

        if len(all_altloc_hier) > 2:
            raise Warning("Only two hierarchies can be compared currently. "
                          "Assuming that if a larger group, that they all overlap if the first two states overlap")

        if hierarchy_all_overlap(all_altloc_hier[0], all_altloc_hier[1], rmsd_distance=0.05):
            print("Altlocs {} are Coincident hierarchies".format(altloc_group))
            coincident_altloc_groups.append(altloc_group)

    return coincident_altloc_groups

def run(pdb,params):

    occupancy_groups = get_occupancy_groups(params)
    altlocs = get_altlocs_from_occupancy_groups(occupancy_groups)
    altloc_groups_with_shared_residues = find_altlocs_with_shared_residues(altlocs, occupancy_groups)

    if len(altloc_groups_with_shared_residues) > 2:
        raise Warning("Current code can only handle two states")

    coincident_altloc_groups = get_coincident_altloc_groups(altloc_groups_with_shared_residues,
                                                            occupancy_groups, pdb)
    return coincident_altloc_groups


# def get_cartesian_grid_points_near_chains(params, inputs, hier):
#     # TODO try replacement with extract xyz
#     xrs_chains = hier.extract_xray_structure(crystal_symmetry=inputs.crystal_symmetry)
#     # print(hier.overall_counts().as_str())
#     sites_cart_chains = xrs_chains.sites_cart()
#
#     # print(params)
#
#     # Calculate the extent of the grid
#     # TODO Test different grid expansion sizes
#     grid_min = flex.double([s - params.options.buffer for s in sites_cart_chains.min()])
#     grid_max = flex.double([s + params.options.buffer for s in sites_cart_chains.max()])
#
#     # TODO Remove dependence on giant.grid
#     grid_near_lig = grid.Grid(grid_spacing=params.options.grid_spacing,
#                               origin=tuple(grid_min),
#                               approx_max=tuple(grid_max))
#
#     return grid_near_lig.cart_points()

        # Genereate hierarcy given altloc and residues


    # for a,b in itertools.combinations(altloc_residue_sets,max_number_altlocs):
    #     if not a.symmetric_difference(b):
    #         print(a, b)





        # For each residue get altlocs of that residue
            # TODO Check whether each group of altlocs are within x RMSD, for all residues in that group
            # TODO Return coinicdent groups of altlocs

    # TODO Count number of altlocs in coincident group (feed to occupancy setting code)

    # Loop over coindident groups
        # Loop over residues in coincident group
            # TODO Get grid cartesian points within buffer of residues
                # TODO Append to overall grid point list for that group
                # TODO Return grids and boolean selections on hiearchies to be used in exhaustive_search.


    # Compare RMSD of altlocs of occupancy groups

    # return ground overlapping altlocs
    # Use length to determine how to dived occupancy

    # return bound overlapping altlocs


if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)