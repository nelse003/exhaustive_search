import os
import numpy as np
import csv

from exhaustive import atom_points_from_sel_string, convex_hull_grid_points, \
    convex_hull_per_residue, convex_hull_from_states
from utils_ccp4 import process_validation_csvs
from utils import u_iso_to_b_fac
from plot import plot_protein_and_selection
from exhaustive import process_refined_pdb_bound_ground_states
from validation import run as validate


def repeat_validate(params):

    # Issue: Can't do selection of lig atoms for each case, can we do it by LIG naming

    # Ligand grid (by convex hull of ligand atoms)
    params.exhaustive.options.ligand_grid_points = True
    print(params.output.out_dir)
    params.output.out_dir = os.path.join(params.output.out_dir, "lig_grid")
    validate(params)
    summary_validation(params)

    # Add plotting of residue selection
    atom_points = atom_points_from_sel_string(pdb=params.input.pdb,
                                              selection_string=
                                              params.exhaustive.options.atom_points_sel_string)

    lig_grid_points = convex_hull_grid_points(atom_points,params)

    plot_protein_and_selection(pdb=params.input.pdb,
                               atom_points=lig_grid_points,
                               plot_filename=os.path.join(
                                   params.output.out_dir,
                                   "lig_grid_points.png"),
                               params=params)

    # Reset
    params.exhaustive.options.ligand_grid_points = False
    params.output.out_dir = os.path.split(params.output.out_dir)[0]

    # Per residue selection of atoms
    params.exhaustive.options.per_residue = True
    params.output.out_dir = os.path.join(params.output.out_dir, "per_residue")
    if os.path.exists(params.output.out_dir):
        for item in os.listdir(params.output.out_dir):
            if item.endswith(".mtz"):
                if not item.startswith("refine"):
                    os.remove(os.path.join(params.output.out_dir, item))
    validate(params)
    summary_validation(params)

    # Add plotting of residue selection
    bound_states, \
    ground_states = process_refined_pdb_bound_ground_states(pdb=params.input.pdb,
                                                            params=params)

    per_residue_points = convex_hull_per_residue(pdb=params.input.pdb,
                                          bound_states=bound_states,
                                          ground_states=ground_states,
                                          params=params)

    plot_protein_and_selection(pdb=params.input.pdb,
                               atom_points=per_residue_points,
                               plot_filename=os.path.join(
                                   params.output.out_dir,
                                   "per_residue_points.png"),
                               params=params)

    #Reset
    params.exhaustive.options.per_residue = False
    params.output.out_dir = os.path.split(params.output.out_dir)[0]
    print(params.output.out_dir)

    # Buffer range (convex hull around occupancy group)
    for buffer_dist in np.arange(0, 1, 0.5):
        params.exhaustive.options.convex_hull_ignore_nearest = False
        params.output.out_dir = os.path.join(params.output.out_dir,
                                             "convex_hull_buffer_{}".format(str(buffer_dist).replace(".", "_")))
        params.exhaustive.options.convex_hull = True
        params.exhaustive.options.buffer = buffer_dist

        if os.path.exists(params.output.out_dir):

            for item in os.listdir(params.output.out_dir):
                if item.endswith(".mtz"):
                    if not item.startswith("refine"):
                        os.remove(os.path.join(params.output.out_dir, item))

        validate(params)
        summary_validation(params)

        buffered_points = convex_hull_from_states(pdb=params.input.pdb,
                                              bound_states=bound_states,
                                              ground_states=ground_states,
                                              params=params)

        plot_protein_and_selection(pdb=params.input.pdb,
                               atom_points=buffered_points,
                               plot_filename=os.path.join(
                                   params.output.out_dir,
                                   "buffer_{}_points.png".format(
                                       str(params.validate.options.set_b).replace('.','_'))),
                               params=params)

        params.output.out_dir = os.path.split(params.output.out_dir)[0]

def summary_validation(params):

    """ Generate a summary csv with validation output

    Output:

    Distance normalised between aimed occupancy and u_iso (b factor)
    - euclidean
    Non normalised distance in B factor and Occ
    """

    min_fofcs, min_occs, min_b_facs, fofcs, occs, b_facs = \
        process_validation_csvs(params.validate.options.start_simul_occ,
                     params.validate.options.end_simul_occ,
                     step=params.validate.options.step_simulation,
                     set_b=params.validate.options.set_b,
                     out_dir=params.output.out_dir,
                     params=params)

    occ_delta = np.abs(np.array(min_occs) - np.array(occs))
    b_delta = np.abs(np.array(min_b_facs) - np.array(b_facs))

    normalised_min_b_fac = (np.array(min_b_facs)
                            - u_iso_to_b_fac(
        params.exhaustive.options.lower_u_iso))/(u_iso_to_b_fac(
        params.exhaustive.options.upper_u_iso)-u_iso_to_b_fac(
        params.exhaustive.options.lower_u_iso))

    normalised_b_fac = (np.array(b_facs)
                            - u_iso_to_b_fac(
        params.exhaustive.options.lower_u_iso))/(u_iso_to_b_fac(
        params.exhaustive.options.upper_u_iso)-u_iso_to_b_fac(
        params.exhaustive.options.lower_u_iso))

    norm_b_delta = np.abs(normalised_min_b_fac - normalised_b_fac)

    mean_occ_delta =  np.mean(occ_delta)
    mean_b_delta =  np.mean(b_delta)
    dst = np.sqrt(occ_delta**2 + norm_b_delta**2 )
    mean_dst = np.mean(dst)

    print(os.path.join(params.output.out_dir, "validation_summary.csv"))

    with open(os.path.join(params.output.out_dir,
                           "validation_summary.csv"),
              'wb') as validation_csv:

        validation_writer = csv.writer(validation_csv, delimiter=',')
        validation_writer.writerow(["mean_occ_delta", "mean_b_delta", "mean_dst"])
        validation_writer.writerow([mean_occ_delta,mean_b_delta,mean_dst])
