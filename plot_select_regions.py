import os
import numpy as np
import matplotlib.pyplot as plt
from exhaustive.exhaustive.utils.convex_hull import convex_hull_per_residue

from exhaustive.exhaustive.utils.select_atoms \
    import process_refined_pdb_bound_ground_states, \
    get_occupancy_group_grid_points

from exhaustive.exhaustive.utils.convex_hull import \
    convex_hull_from_states, atom_points_from_sel_string, \
    convex_hull_grid_points

from exhaustive.exhaustive.utils.utils_ccp4 import expand_array


# params =  master_phil.extract()
# params.input.xtal_name = "NUDT22A-x0182"
# params.input.in_path = os.path.join("/hdlocal/home/enelson/exhaustive_search_data/validation_NUDT22/FMOPL000622a_DSPL/NUDT22A-x0182")
# params.input.mtz = os.path.join(params.input.in_path, "refine.mtz")
# params.input.pdb = os.path.join(params.input.in_path,"refine.pdb")
# params.output.out_dir = os.path.realpath("/hdlocal/home/enelson/exhaustive_search_data/validation_NUDT22/FMOPL000622a_DSPL/NUDT22A-x0182")
# params.exhaustive.options.atom_points_sel_string = "(chain B and altid C and resid 1) or (chain B and altid D resid 1)"

def plot_protein_region(params, lig_grid=False, all=True, per_residue=True,
                        convex_hull=False, boxes=False, lig_atoms=False):

    bound_states, \
    ground_states = process_refined_pdb_bound_ground_states(pdb=params.input.pdb,
                                                            params=params)

    #atom_points_within_states(pdb=params.input.pdb, bound_states=bound_states, ground_states=ground_states)

    convex_hull_points = convex_hull_from_states(pdb=params.input.pdb,
                                          bound_states=bound_states,
                                          ground_states=ground_states,
                                          params=params)

    convex_hull_array = np.array(list(convex_hull_points))

    lig_atom_points = atom_points_from_sel_string(pdb=params.input.pdb,
                                              selection_string=
                                              params.exhaustive.options.atom_points_sel_string)

    lig_atom_array = np.array(list(lig_atom_points))

    lig_grid_points = convex_hull_grid_points(lig_atom_points, params=params)

    lig_grid_point_array = np.array(list(lig_grid_points))

    box_points = get_occupancy_group_grid_points(pdb=params.input.pdb,
                                                  bound_states=bound_states,
                                                  ground_states=ground_states,
                                                  params=params)
    box_array = np.array(list(box_points))

    per_residue_convex_hull_points = convex_hull_per_residue(pdb=params.input.pdb,
                            bound_states=bound_states,
                            ground_states=ground_states,
                            params=params)

    per_residue_convex_hulls_array = np.array(list(per_residue_convex_hull_points))

    #whole protein?

    all_atoms =atom_points_from_sel_string(pdb=params.input.pdb, selection_string='all')
    # all_points = convex_hull_grid_points(all_atoms, params=params)
    # print(all_points.size())

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.axis('equal')

    if lig_grid:
        # # Ligand grid points
        x,y,z = expand_array(lig_grid_point_array)
        ax.scatter(x,y,z, marker='.', color='k')

    if per_residue:
        # Convex Hull per residue
        x,y,z = expand_array(per_residue_convex_hulls_array)
        ax.scatter(x,y,z, marker='.', color='m')

    if convex_hull:
        # Convex_hull
        x,y,z = expand_array(convex_hull_array)
        ax.scatter(x,y,z, marker='.', color='g')

    if boxes:
        # # Boxes around lig of interest
        x,y,z = expand_array(box_array)
        ax.scatter(x,y,z, marker='.', color='c')

    if all:
        # All atoms
        x,y,z = expand_array(np.array(list(all_atoms)))
        ax.scatter(x,y,z, marker='.', color='y')

    if lig_atoms:
        # # Ligand atom points
        x,y,z = expand_array(lig_atom_array)
        ax.scatter(x,y,z, marker='o', color='r')

    plt.savefig(filename=os.path.join(params.output.out_dir,"region_selected.png"), dpi=300)
