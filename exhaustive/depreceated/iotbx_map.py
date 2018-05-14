
# This can be used to write ccp4 maps directly. Was used in exhaustive.py,
# where inputs have been read in, and fofc-map generated

    iotbx.ccp4_map.write_ccp4_map(
        file_name ="testing_iotbx_{}_{}.ccp4".format(bound_occupancy, u_iso),
        unit_cell = inputs.crystal_symmetry.unit_cell(),
        space_group = inputs.crystal_symmetry.space_group()
        map_data = fofc_map
    )

# For converting fft_map to ccp4 map. Works on map, not ofn coefficents so may onyl show are used?

fft_map.as_ccp4_map(file_name="testing_{}_{}.ccp4".format(bound_occupancy, u_iso))