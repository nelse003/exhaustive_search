import iotbx.pdb
import numpy as np

from scipy.spatial import ConvexHull
from scipy.spatial import distance

from scitbx.array_family import flex

from select_atoms import get_occupancy_groups

def convex_hull_from_states(pdb, bound_states, ground_states, params):

    """ Convex hull selection of bound & ground states 
    
    Given bound & ground states and pdb file, determine a grid of points which
    are defined by a convex hull around the atoms involved in the bound and 
    ground states. 
    
    Allow extension of the convex hull by a buffer distance, 
    supplied via params.exhaustive.options.buffer. The buffer extends 
    in the direction of the nearest atom not in the bound or ground states.
    By default (params.exhaustive.options.convex_hull_ignore_nearest = False)
    the buffer is at most 1/2 the distance to the nearest atom. When 
    overridden the buffer will continue in the direction of the nearest atom 
    up to the supplied buffer distance.
    
    :return: flex.vec3_double containing the x,y,z coordinates of points
              within the convex hull
      
    """

    # Get atom points that come from states
    atom_points = atom_points_within_states(pdb, bound_states, ground_states)

    # If no buffer
    if params.exhaustive.options.buffer == 0:
        return convex_hull_grid_points(atom_points, params)

    # Select atoms that are involved in the occupancy groups
    # and remove these from the distance matrix
    # This will mean that atoms selected as close by are
    # not already in that state
    atoms_not_in_occ_group = atoms_not_in_occ_groups(pdb, params)
    atoms_not_in_occ_group_xyz = atoms_not_in_occ_group.extract_xyz()

    # generate a convex hull
    hull = ConvexHull(atom_points)

    # Find atoms closest to the convex hull, which are not part of the convex hull
    buffer_points = []
    for vertex in hull.vertices:

        print(list(atom_points)[vertex])
        vertex_atom_point = [list(atom_points)[vertex]]
        dist_matrix = distance.cdist(vertex_atom_point, atoms_not_in_occ_group_xyz)
        min_dist = dist_matrix.min()
        index_atom = dist_matrix.argmin()

        if min_dist == 0:
            raise ValueError("Minimal distance to nearest atom is 0.0")

        nearest_atom_point = np.array(atoms_not_in_occ_group_xyz[index_atom])

        print(nearest_atom_point)
        print(type(nearest_atom_point))
        print(vertex_atom_point)
        print(type(vertex_atom_point))

        if not params.exhaustive.options.convex_hull_ignore_nearest \
            and min_dist/2 >= params.exhaustive.options.buffer:

            buffer_point = vertex_atom_point + \
                           (params.exhaustive.options.buffer/min_dist)*(nearest_atom_point - vertex_atom_point)

            print("Buffer point at buffer distance")

        elif not params.exhaustive.options.convex_hull_ignore_nearest \
            and min_dist/2 < params.exhaustive.options.buffer:

            buffer_point = vertex_atom_point + 0.5*(nearest_atom_point - vertex_atom_point)
            print("Buffer point at mid point to nearest atom")

        else:
            buffer_point = vertex_atom_point + \
                           params.exhaustive.options.buffer * \
                           (nearest_atom_point - vertex_atom_point)/min_dist
            print("Buffer point ignores nearest atom")

        buffer_points.append(buffer_point)

    buffer_points = np.concatenate(buffer_points)

    return convex_hull_grid_points(buffer_points, params)

def convex_hull_grid_points(atom_points,params):

    """
    Given points define grid points in convex hull
    
    Given points, either nd.array or flex.vec3_double, 
    return a grid of points within a convex hull 
    defined by those points 
    
    :param atom_points: 
    :param params: 
    :return: flex.vec3_double 
    """

    hull = ConvexHull(atom_points)

    if  isinstance(atom_points,flex.vec3_double):
        grid_min = flex.double(atom_points.min())
        grid_max = flex.double(atom_points.max())
    elif isinstance(atom_points,np.ndarray):
        grid_min = flex.double(atom_points.min(axis=0))
        grid_max = flex.double(atom_points.max(axis=0))
    else:
        raise TypeError("Convex hull grids need either numpy array, or flex.vec3_double")

    grid_from_selection = grid.Grid(
        grid_spacing=params.exhaustive.options.grid_spacing,
        origin=tuple(grid_min),
        approx_max=tuple(grid_max))

    grid_points = grid_to_np_array(grid_from_selection)

    points_in_hull = []

    for point in grid_points:
        if point_in_hull(point, hull):
            points_in_hull.append(point)

    print("Points in hull: {}".format(len(points_in_hull)))

    return flex.vec3_double(list(np.stack(points_in_hull)))


def atom_points_within_states(pdb, bound_states, ground_states):

    """ Given states return xyz of atoms involved
    
    
    :param pdb: 
    :param bound_states: 
    :param ground_states: 
    :return: flex.vec3_double
    """

    states = bound_states + ground_states
    pdb_in = iotbx.pdb.hierarchy.input(pdb)
    pdb_atoms = pdb_in.hierarchy.atoms()
    atom_points = flex.vec3_double()
    all_selected_atoms = []

    for state in states:
        selection = state[0]
        selected_atoms = pdb_atoms.select(selection)
        all_selected_atoms.append(selected_atoms)
        sites_cart = selected_atoms.extract_xyz()
        atom_points = atom_points.concatenate(sites_cart)

    return atom_points

def atoms_not_in_occ_groups(pdb, params):

    """Select atoms that are not involved in the occupancy groups
    
    :returns: iotbx.pdb.hierarchy.af_shared_atom
    
    """

    occupancy_groups = get_occupancy_groups(pdb, params)
    selection_string_list = []
    for group in occupancy_groups:
        for item in group:
            for residue in item:
                selection_string = "(chain {} " \
                                   "and altid {} " \
                                   "and resid {})".format(residue['chain'],
                                                          residue['altloc'],
                                                          residue['resseq'])
                selection_string_list.append(selection_string)

    selection_string = " or ".join(selection_string_list)
    not_selection_string ="not ({})".format(selection_string)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()
    remove_atoms_sel = sel_cache.selection(not_selection_string)
    removed_hier = pdb_in.hierarchy.select(remove_atoms_sel)

    atoms_not_in_occ_group = removed_hier.atoms()
    return removed_hier.atoms()


def point_in_hull(point, hull, tolerance=1e-12):
    """
    https://stackoverflow.com/questions/16750618/
    whats-an-efficient-way-to-find-if-a-
    point-lies-in-the-convex-hull-of-a-point-cl
    """

    return all(
        (np.dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)

def grid_to_np_array(grid):

    array = np.zeros((grid.cart_points().all()[0],3))
    for i, point in enumerate(grid.cart_points()):
        array[i] = point

    return array