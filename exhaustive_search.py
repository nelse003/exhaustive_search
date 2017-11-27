#########################################################################
#  Imports
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

import giant.xray.edstats as ed
########################################################################

def compute_maps(fmodel, crystal_gridding, map_type):
    map_coefficients = map_tools.electron_density_map(
        fmodel = fmodel).map_coefficients(
        map_type         = map_type,
        isotropize       = True,
        fill_missing     = False)
    fft_map = cctbx.miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = map_coefficients)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded(), map_coefficients

#########################################################################
#
#   Cmd Line arguments to run the program:
#   script filepath_to_pdb filepath_to_mtz
#
#########################################################################

# Loop over atom position , b_factor and occupancy

def loop_over_atoms(sites_frac, sel_lig, xrs, atoms, fmodel, crystal_gridding):

    # TODO Swap to only generate fraction sites related to selected ligand (i.e remove if statement)
    # TODO Profile the speed of the loop

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

                        ## Write out MTZ
                        # mtz_dataset = fofc.as_mtz_dataset(column_root_label="FOFCWT")
                        # mtz_object = mtz_dataset.mtz_object()
                        # mtz_object.write(file_name = "{}_{}.mtz".format(prefix,occupancy))
                        ##

                        fofc_value = fofc_map.eight_point_interpolation(site_frac)
                        print(occupancy, name, "%8.3f" % (fofc_value))
                        row = [xrs_dc.scatterers()[i].occupancy, xrs_dc.scatterers()[i].u_iso, fofc_value]
                        writer.writerow(row)
                        sys.stdout.flush()

# TODO Turn the loop statements into generator, reduce repeating code
def loop_over_residues_edstats(sites_frac, sel_lig, xrs, atoms, fmodel, crystal_gridding):

    # TODO get residue identifier?
    with open('LIG.csv'.format(),'w') as f1:
        writer = csv.writer(f1, delimiter=',', lineterminator='\n')
        # currently loop over rough occupancy range for initial testing
        for occupancy in numpy.arange(0, 1.01, 0.05):
            for u_iso in numpy.arange(0.25, 1.2, 0.05):
                xrs_dc = xrs.deep_copy_scatterers()

                # TODO Set occupancy using whole residue group rather than for loop
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



def loop_over_atoms_find_fofc_at_multiple_sites(sites_frac, sel_lig, xrs, atoms, fmodel, crystal_gridding):

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

# TODO Get mean B factor of crystal: Use as input to loop over occupancy with fixed B-Factor
# Use edstats?

# TODO Get B factor of surrounding residues: Use as input to loop over occupancy with fixed B-Factor
# Use iotbx heirarchy?

def cmd_run(args, out=sys.stdout):

    # Read input PDB and reflection files
    inputs = mmtbx.utils.process_command_line_args(args = args)
    reflection_files = inputs.reflection_files
    rfs = reflection_file_utils.reflection_file_server(
        crystal_symmetry = inputs.crystal_symmetry,
        force_symmetry   = True,
        reflection_files = reflection_files,
        err              = StringIO())
    determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
        reflection_file_server = rfs,
        keep_going             = True,
        log                    = StringIO())
    pdb_inp = iotbx.pdb.input(file_name = inputs.pdb_file_names[0])
    ph = pdb_inp.construct_hierarchy()
    asc = ph.atom_selection_cache()

    # Select ligand atoms
    sel_lig = asc.selection("resname LIG")
    xrs = ph.extract_xray_structure(
        crystal_symmetry = inputs.crystal_symmetry)
    # Extract Fobs and free-r flags
    f_obs = determined_data_and_flags.f_obs
    r_free_flags = determined_data_and_flags.r_free_flags
    # Define map grididng
    crystal_gridding = f_obs.crystal_gridding(
        d_min             = f_obs.d_min(),
        symmetry_flags    = maptbx.use_space_group_symmetry,
        resolution_factor = 1./4)
    # Define fmodel
    mask_params = mmtbx.masks.mask_master_params.extract()
    mask_params.ignore_hydrogens=False
    mask_params.ignore_zero_occupancy_atoms=False
    fmodel = mmtbx.f_model.manager(
        f_obs          = f_obs,
        r_free_flags   = r_free_flags,
        mask_params    = mask_params,
        xray_structure = xrs)
    fmodel.update_all_scales()
    print("r_work: {0} r_free: {1}".format(fmodel.r_work(), fmodel.r_free()))
    # Show map values at ligand atom centers
    sites_frac = xrs.sites_frac()
    atoms = list(ph.atoms())
    sel = xrs.hd_selection()

    # Generate results in an output directory

    output_folder = "output_a"
    output_path = os.path.join(os.getcwd(),output_folder)
    print(output_path)

    if not os.path.exists(output_path):
        os.mkdir(output_path)
    # TODO Swap out the change directory to write out for each function writing out.
    os.chdir(output_folder)

    #loop_over_atoms_find_fofc_at_multiple_sites(sites_frac, sel_lig, xrs, atoms, fmodel, crystal_gridding)
    loop_over_residues_edstats(sites_frac, sel_lig, xrs, atoms, fmodel, crystal_gridding)
    #loop_over_atoms(sites_frac, sel_lig, xrs, atoms, fmodel, crystal_gridding)

if(__name__ == "__main__"):
    cmd_run(args = sys.argv[1:])