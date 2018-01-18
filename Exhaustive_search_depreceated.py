from __future__ import print_function
from __future__ import division
import sys
from cStringIO import StringIO
import mmtbx.f_model
import mmtbx.utils
from mmtbx import map_tools
from iotbx import reflection_file_utils
import iotbx.pdb
from cctbx import maptbx
import cctbx.miller
import mmtbx.masks
import numpy
import os
import csv
from scitbx.array_family import flex
import giant.xray.edstats as ed
import giant.grid as grid
import numpy as np
import matplotlib.pyplot as plt

# Loop over atom position , b_factor and occupancy

def loop_over_atoms(sites_frac, sel_lig, xrs, atoms, fmodel, crystal_gridding):

    for i, site_frac in enumerate(sites_frac):
        if (sel_lig[i]):
            with open('{}.csv'.format(atoms[i].format_atom_record()[:11].replace(" ", "_")), 'w') as f1:
                writer = csv.writer(f1, delimiter=',', lineterminator='\n')
                # currently loop over rough occupancy range for initial testing
                for occupancy in numpy.arange(0, 1.01, 0.05):
                    for u_iso in numpy.arange(0.25, 1.2, 0.05):
                        xrs_dc = xrs.deep_copy_scatterers()
                        # vary occupancy of scatterer?
                        xrs_dc.scatterers()[i].occupancy = occupancy
                        xrs_dc.scatterers()[i].u_iso = u_iso
                        print(xrs_dc.scatterers()[i].u_iso, xrs_dc.scatterers()[i].occupancy)
                        fmodel.update_xray_structure(
                            xray_structure=xrs_dc,
                            update_f_calc=True)
                        fofc_map, fofc = compute_maps(
                            fmodel=fmodel,
                            crystal_gridding=crystal_gridding,
                            map_type="mFo-DFc")
                        name = atoms[i].format_atom_record()[:28]
                        prefix = "_".join(name.split())

                        fofc_value = fofc_map.eight_point_interpolation(site_frac)
                        print(occupancy, name, "%8.3f" % (fofc_value))
                        row = [xrs_dc.scatterers()[i].occupancy, xrs_dc.scatterers()[i].u_iso, fofc_value]
                        writer.writerow(row)
                        sys.stdout.flush()

def loop_over_residues_edstats(sites_frac, sel_lig, xrs, atoms, fmodel, crystal_gridding):

    with open('LIG.csv'.format(),'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        # currently loop over rough occupancy range for initial testing
        for occupancy in numpy.arange(0, 1.01, 0.05):
            for u_iso in numpy.arange(0.25, 1.2, 0.05):
                xrs_dc = xrs.deep_copy_scatterers()

                # Change Occupancy and B factor for all atoms in selected ligand at the same time
                for i, site_frac in enumerate(sites_frac):
                    if (sel_lig[i]):
                        xrs_dc.scatterers()[i].occupancy = occupancy
                        xrs_dc.scatterers()[i].u_iso = u_iso
                        print(xrs_dc.scatterers()[i].u_iso, xrs_dc.scatterers()[i].occupancy)

                fmodel.update_xray_structure(
                    xray_structure=xrs_dc,
                    update_f_calc=True)
                """fofc_map, fofc = compute_maps(
                    fmodel=fmodel,
                    crystal_gridding=crystal_gridding,
                    map_type="mFo-DFc")
                name = atoms[i].format_atom_record()[:28]
                fofc_value = fofc_map.eight_point_interpolation(site_frac)
                print(occupancy, name, "%8.3f" % (fofc_value))"""

                # TODO Generate 2fofc maps for Edstats?

                two_fofc_map, two_fofc = compute_maps(
                    fmodel=fmodel,
                    crystal_gridding=crystal_gridding,
                    map_type="2mFo-DFc")

                print(two_fofc_map)
                print(type(two_fofc_map))
                print(two_fofc)
                print(type(two_fofc))
                # TODO use edstats directely on outputted map
                # TODO Or replace edstats with just RSR generation function/ program?

                mtz_dataset = two_fofc.as_mtz_dataset(column_root_label="2FOFCWT")
                mtz_object = mtz_dataset.mtz_object()
                mtz_object.write(file_name = "LIG_{}.mtz".format(occupancy))

                # TODO Generalise map loaded in?
                # Calclulate RSR value?
                print(os.getcwd())
                edstats, summary = ed.score_file_with_edstats("LIG_{}.mtz".format(occupancy),
                                                              "/hdlocal/home/enelson/Dropbox/DPhil/exhaustive_search/refine_1.pdb")

                print(edstats.scores)

                # TODO Utilise edstat function to select residue group for ligand (Check the way nick uses this)

                # RSR score for ligand
                RSR_LIG = edstats.scores.loc['Ra', (slice(None), 'F', slice(None), slice(None))]
                print(RSR_LIG)

                row = [xrs_dc.scatterers()[i].occupancy, xrs_dc.scatterers()[i].u_iso, fofc_value]
                writer.writerow(row)
                sys.stdout.flush()

def loop_over_residues_sum_fofc(sites_frac, sel_lig, xrs, xrs_lig, atoms, fmodel, crystal_gridding):

    sites_cart_lig = xrs_lig.sites_cart()

    # Calculate the extent of the grid
    buffer = 5
    grid_min = flex.double([s-buffer for s in sites_cart_lig.min()])
    grid_max = flex.double([s+buffer for s in sites_cart_lig.max()])

    grid_near_lig = grid.Grid(grid_spacing = 0.5,
                              origin = tuple(grid_min),
                              approx_max = tuple(grid_max))

    sites_cart_near_lig = grid_near_lig.cart_points()
    print(len(list(sites_cart_lig)))
    print(len(list(sites_cart_near_lig)))

    with open('mean_5_buffer_bound_only_no_covary.csv'.format(),'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        # currently loop over rough occupancy range for initial testing
        for occupancy in numpy.arange(0, 1.01, 0.05):
            for u_iso in numpy.arange(0.25, 1.2, 0.05):
                xrs_dc = xrs.deep_copy_scatterers()

                # Change Occupancy and B factor for all atoms in selected ligand at the same time
                for i, site_frac in enumerate(sites_frac):
                    if (sel_lig[i]):
                        xrs_dc.scatterers()[i].occupancy = occupancy
                        xrs_dc.scatterers()[i].u_iso = u_iso
                        #print(xrs_dc.scatterers()[i].u_iso, xrs_dc.scatterers()[i].occupancy)

                fmodel.update_xray_structure(
                    xray_structure=xrs_dc,
                    update_f_calc=True)
                fofc_map, fofc = compute_maps(
                    fmodel=fmodel,
                    crystal_gridding=crystal_gridding,
                    map_type="mFo-DFc")
                name = atoms[i].format_atom_record()[:28]

                sum_abs_fofc_value = 0
                for i, site_cart_near_lig in enumerate(sites_cart_near_lig):
                    site_frac_near_lig = inputs.crystal_symmetry.unit_cell().fractionalize(site_cart_near_lig)
                    fofc_value = fofc_map.eight_point_interpolation(site_frac_near_lig)
                    sum_abs_fofc_value += abs(fofc_value)

                print(occupancy, "%8.3f" % (sum_abs_fofc_value))

                mean_abs_fofc_value = sum_abs_fofc_value/len(list(sites_cart_near_lig))

                row = [occupancy, u_iso, mean_abs_fofc_value]
                writer.writerow(row)
                sys.stdout.flush()

# TODO Check fractionalisation code in this function!!
def get_mean_fofc(xrs, sites_frac, fmodel, crystal_gridding):

    sites_cart = xrs.sites_cart()

    # Get grid for whole ligand
    buffer = 5
    grid_min = flex.double([s-buffer for s in sites_cart.min()])
    grid_max = flex.double([s+buffer for s in sites_cart.max()])

    # TODO Remove dependence on giant.grid
    grid_protein = grid.Grid(grid_spacing = 0.25,
                              origin = tuple(grid_min),
                              approx_max = tuple(grid_max))

    grid_sites_cart = grid_protein.cart_points()
    grid_sites_frac = grid_protein.frac_points()

    fofc_map, fofc = compute_maps(
        fmodel=fmodel,
        crystal_gridding=crystal_gridding,
        map_type="mFo-DFc")

    sum_fofc_value = 0
    sum_abs_fofc_value = 0

    for i, site_frac in enumerate(sites_frac):
        fofc_value = fofc_map.eight_point_interpolation(site_frac)
        sum_fofc_value += fofc_value
        sum_abs_fofc_value += abs(fofc_value)

    # TODO Find a better way to get this over grid mean fo-fc
    sum_grid_fofc = 0
    sum_grid_abs_fofc = 0
    for site in grid_sites_frac:
        fofc_value = fofc_map.eight_point_interpolation(site)
        sum_grid_fofc += fofc_value
        sum_grid_abs_fofc += abs(fofc_value)


    print("Mean Fo-Fc {}".format(sum_fofc_value / len(sites_frac)))
    print("Mean |Fo-Fc| {}".format(sum_abs_fofc_value / len(sites_frac)))

    print("Mean grid Fo-Fc {}".format(sum_grid_fofc / len(grid_sites)))
    print("Mean grid |Fo-Fc| {}".format(sum_grid_abs_fofc / len(grid_sites)))

    return sum_fofc_value / len(sites_frac), \
           sum_abs_fofc_value / len(sites_frac), \
           sum_grid_fofc/len(grid_sites), \
           sum_grid_abs_fofc/len(grid_sites)

def loop_over_atoms_find_fofc_at_multiple_sites(sites_frac, sel_lig, xrs, xrs_lig, atoms, fmodel, crystal_gridding):

    #TODO Select a range of sites to loop over to sample the map calculated for each position

    for i, site_frac in enumerate(sites_frac):
        with open('occ_050_uiso_035_.csv', 'w') as f1:
            writer = csv.writer(f1, delimiter=',', lineterminator='\n')
            # currently loop over rough occupancy range for initial testing
            #for occupancy in numpy.arange(0, 1.01, 0.05):
            #    for u_iso in numpy.arange(0.25, 1.2, 0.05):
            xrs_dc = xrs.deep_copy_scatterers()
            # vary occupancy of scatterer?
            xrs_dc.scatterers()[i].occupancy = 0.5
            xrs_dc.scatterers()[i].u_iso = 0.35
            fmodel.update_xray_structure(
                xray_structure=xrs_dc,
                update_f_calc=True)
            fofc_map, fofc = compute_maps(
                fmodel=fmodel,
                crystal_gridding=crystal_gridding,
                map_type="mFo-DFc")
            name = atoms[i].format_atom_record()[:28]
            print(name,xrs_dc.scatterers()[i].u_iso, xrs_dc.scatterers()[i].occupancy)
            all_sites_fofc = [name]
            for j, site_frac_X in enumerate(sites_frac):
                fofc_value = fofc_map.eight_point_interpolation(site_frac_X)
                all_sites_fofc.append(fofc_value)
                writer.writerow(all_sites_fofc)
            sys.stdout.flush()


def fofc_between_two_atoms(site_cart_atom_1, site_cart_atom_2, fmodel,
                           crystal_gridding, inputs, buffer=2, spacing=0.01):
    """ |Fo-Fc| value plotted for points between two atoms, extended by a buffer"""

    site_cart_atom_1_np = np.asarray(site_cart_atom_1)
    site_cart_atom_2_np = np.asarray(site_cart_atom_2)

    length = abs(np.sqrt(site_cart_atom_2_np.dot(site_cart_atom_2_np)) - np.sqrt(site_cart_atom_1_np.dot(site_cart_atom_1_np)))
    intervals = length/spacing
    sampling = np.floor(intervals)
    sample_vector = (site_cart_atom_2_np-site_cart_atom_1_np)/sampling
    adjusted_sampling =length/sampling

    buffer_intervals = int(buffer/adjusted_sampling)

    fofc_map, fofc = compute_maps(
        fmodel=fmodel,
        crystal_gridding=crystal_gridding,
        map_type="mFo-DFc")

    all_fofc=[]
    line_coords = adjusted_sampling * np.asarray(range (0 - buffer_intervals, int(sampling) + buffer_intervals))

    for i in range(0 - buffer_intervals, int(sampling) + buffer_intervals):
        location = sample_vector * i + site_cart_atom_1
        location_frac = inputs.crystal_symmetry.unit_cell().fractionalize(location)
        fofc_value = fofc_map.eight_point_interpolation(location_frac)
        all_fofc.append(fofc_value)

    return line_coords, all_fofc

def fofc_line_plot(line_coords, all_fofc, plot_name, buffer=2):

    plt.plot(line_coords,all_fofc)
    plt.xlabel(" Distance from Atom 1")
    plt.ylabel(" Fo - Fc value ")
    plt.savefig(plot_name)
    plt.close()

# TODO Get mean B factor of crystal: Use as input to loop over occupancy with fixed B-Factor
# Use edstats?

# TODO Get B factor of surrounding residues: Use as input to loop over occupancy with fixed B-Factor
# Use iotbx heirarchy?

# From exhaustive_search run()

    ## Fo-Fc Line plotting

    # fofc_map, fofc = compute_maps(
    #     fmodel=fmodel,
    #     crystal_gridding=crystal_gridding,
    #     map_type="mFo-DFc")
    #
    # lig_cart = xrs_lig.sites_cart()
    # lig_frac = xrs_lig.sites_frac()
    #
    # print(type(lig_cart[0]))
    # print(fofc_map.eight_point_interpolation(lig_cart[0]), fofc_map.eight_point_interpolation(lig_frac[0]))
    # print(fofc_map.eight_point_interpolation(inputs.crystal_symmetry.unit_cell().fractionalize(lig_cart[0])))
    #
    # for i in range(1,len(lig_cart)):
    #     line_coords, all_fofc = fofc_between_two_atoms(lig_cart[i-1], lig_cart[i], fmodel, crystal_gridding, inputs, buffer=2)
    #     plot_name = "Atom_{}_to_{}".format(i-1,i)
    #     fofc_line_plot(line_coords, all_fofc, plot_name, buffer=2)