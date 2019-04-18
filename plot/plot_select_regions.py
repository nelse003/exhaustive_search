import os

import matplotlib.pyplot as plt
import numpy as np

from exhaustive import atom_points_from_sel_string
from exhaustive import convex_hull_from_states
from exhaustive import convex_hull_grid_points
from exhaustive import convex_hull_per_residue
from exhaustive import get_occupancy_group_grid_points
from exhaustive import process_refined_pdb_bound_ground_states
from utils import expand_array


def plot_protein_region(params, lig_grid=False, all=True, per_residue=True,
                        convex_hull=False, boxes=False, lig_atoms=False):
    bound_states, \
    ground_states = process_refined_pdb_bound_ground_states(pdb=params.input.pdb,
                                                            params=params)

    # atom_points_within_states(pdb=params.input.pdb, bound_states=bound_states, ground_states=ground_states)

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

    # whole protein?

    all_atoms = atom_points_from_sel_string(pdb=params.input.pdb, selection_string='all')
    # all_points = convex_hull_grid_points(all_atoms, params=params)
    # print(all_points.size())

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.axis('equal')

    if lig_grid:
        # # Ligand grid points
        x, y, z = expand_array(lig_grid_point_array)
        ax.scatter(x, y, z, marker='.', color='k')

    if per_residue:
        # Convex Hull per residue
        x, y, z = expand_array(per_residue_convex_hulls_array)
        ax.scatter(x, y, z, marker='.', color='m')

    if convex_hull:
        # Convex_hull
        x, y, z = expand_array(convex_hull_array)
        ax.scatter(x, y, z, marker='.', color='g')

    if boxes:
        # # Boxes around lig of interest
        x, y, z = expand_array(box_array)
        ax.scatter(x, y, z, marker='.', color='c')

    if all:
        # All atoms
        x, y, z = expand_array(np.array(list(all_atoms)))
        ax.scatter(x, y, z, marker='.', color='y')

    if lig_atoms:
        # # Ligand atom points
        x, y, z = expand_array(lig_atom_array)
        ax.scatter(x, y, z, marker='o', color='r')

    plt.savefig(filename=os.path.join(params.output.out_dir, "region_selected.png"), dpi=300)
