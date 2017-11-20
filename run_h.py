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

import os
print os.environ["PYTHONPATH"]

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
  sel_H = asc.selection("protein and element H")
  xrs = ph.extract_xray_structure(
    crystal_symmetry = inputs.crystal_symmetry)
  # Extract Fobs and free-r flags
  f_obs = determined_data_and_flags.f_obs
  r_free_flags = determined_data_and_flags.r_free_flags
  # Define map gridding
  crystal_gridding = f_obs.crystal_gridding(
    d_min             = f_obs.d_min(),
    symmetry_flags = maptbx.use_space_group_symmetry,
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
  print "r_work: %6.4f r_free: %6.4f"%(fmodel.r_work(), fmodel.r_free())
  # Show map values at H atom centers
  sites_frac = xrs.sites_frac()
  atoms = list(ph.atoms())
  # Select all the hydrogen atoms
  sel = xrs.hd_selection()
  for i, site_frac in enumerate(sites_frac):
    if(sel_H[i]):
      #Copy the scatteres from the model
      xrs_dc = xrs.deep_copy_scatterers()
      # Set the occupancy of the each atom to zero?        
      xrs_dc.scatterers()[i].occupancy=0
      # Update the model with scatterers 
      fmodel.update_xray_structure(
        xray_structure = xrs_dc,
        update_f_calc  = True)
      # calculate map and map coefficents
      fofc_map, fofc = compute_maps(
        fmodel           = fmodel, 
        crystal_gridding = crystal_gridding,
        map_type         = "mFo-DFc")
      # get identifying elements in PDB format
      name = atoms[i].format_atom_record()[:28]
      prefix = "_".join(name.split())
      ##
      # Write fofc to mtzfile (One for each H atom) 
      mtz_dataset = fofc.as_mtz_dataset(column_root_label="FOFCWT")
      mtz_object = mtz_dataset.mtz_object()
      mtz_object.write(file_name = "%s.mtz"%prefix)
      ##
      # Interpoalte map to get value at this point
      fofc_value = fofc_map.eight_point_interpolation(site_frac)
      print name, "%8.3f"%(fofc_value)
      sys.stdout.flush()
  

if(__name__ == "__main__"):
  cmd_run(args = sys.argv[1:])

