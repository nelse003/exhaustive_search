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

def get_occupancy_group_grid_points(pdb, params):

    all_occupancy_group_cart_points = get_list_occupancy_group_grid_points(pdb, params)

    occupancy_group_cart_points = flex.vec3_double()
    for occupancy_group_points in all_occupancy_group_cart_points:
        occupancy_group_cart_points = occupancy_group_cart_points.concatenate(occupancy_group_points)

    return occupancy_group_cart_points

def get_list_occupancy_group_grid_points(pdb, params):

    """Return grid points for each occupancy group in a list"""

    all_occupancy_group_cart_points_list = []
    # Get grid points per residue basis, as residues could be seperated in space
    for altloc_hier in get_altloc_hier(pdb):

        print(pdb)

        occupancy_group_cart_points = flex.vec3_double()

        for chain in altloc_hier.only_model().chains():
            for residue_group in chain.residue_groups():

                # TODO Check: Using extract xyz instead of xrs (will this be sufficent?)
                sites_residue_cart = residue_group.atoms().extract_xyz()
                grid_min = flex.double([s - params.options.buffer for s in sites_residue_cart.min()])
                grid_max = flex.double([s + params.options.buffer for s in sites_residue_cart.max()])

                grid_residue = grid.Grid(grid_spacing = params.options.grid_spacing,
                              origin = tuple(grid_min),
                              approx_max = tuple(grid_max))

                occupancy_group_cart_points = occupancy_group_cart_points.concatenate(grid_residue.cart_points())

        print(len(occupancy_group_cart_points))
        all_occupancy_group_cart_points_list.append(occupancy_group_cart_points)

    return all_occupancy_group_cart_points_list

##############################################
#From exhaustive_search run()
##############################################

# Fo-Fc Line plotting

fofc_map, fofc = compute_maps(
    fmodel=fmodel,
    crystal_gridding=crystal_gridding,
    map_type="mFo-DFc")

lig_cart = xrs_lig.sites_cart()
lig_frac = xrs_lig.sites_frac()

print(type(lig_cart[0]))
print(fofc_map.eight_point_interpolation(lig_cart[0]), fofc_map.eight_point_interpolation(lig_frac[0]))
print(fofc_map.eight_point_interpolation(inputs.crystal_symmetry.unit_cell().fractionalize(lig_cart[0])))

for i in range(1,len(lig_cart)):
    line_coords, all_fofc = fofc_between_two_atoms(lig_cart[i-1], lig_cart[i], fmodel, crystal_gridding, inputs, buffer=2)
    plot_name = "Atom_{}_to_{}".format(i-1,i)
    fofc_line_plot(line_coords, all_fofc, plot_name, buffer=2)

#################################################
#Depreceated functions from select_occupancy_groups
#################################################

def get_largest_concident_group(pdb, residue_altloc_dict, params):

    coincident = get_coincident_group(pdb, residue_altloc_dict, params)

    for residue_chain in residue_altloc_dict.keys():
        residue = residue_chain[0]
        chain = residue_chain[1]
        altloc_group = residue_altloc_dict.get(residue_chain)
        non_coincident_altlocs = list(altloc_group)
        for conincident_altloc_group in coincident:
            if conincident_altloc_group[1] == residue and conincident_altloc_group[2] == chain:
                for altloc in conincident_altloc_group[0]:
                    non_coincident_altlocs.remove(altloc)
        if non_coincident_altlocs:
            raise ValueError("More than two groups of altlocs are present. This is currently not handled")

######################################################################
# Selection using altlocs where all residues are shared, and they overlap
######################################################################

# Use occupancy groups to generate selection of atoms associated to bound & ground structure

# Return an atom selection on the protein hier to be used as selection
# Provide whether sufficent for co-varying

# Need to consider generation of grid (for summing FoFc) for seperated components.

###########################################################################
# Used in generation of cartesian point selection
###########################################################################
def generate_altloc_hiearchy(altloc_residue_dict, pdb):

    pdb_in = hierarchy.input(pdb)
    altloc_selection = generate_altloc_selection(altloc_residue_dict, pdb)
    altloc_hierarchy = pdb_in.hierarchy.select(altloc_selection)

    print_hier_atoms(altloc_hierarchy)

    return altloc_hierarchy

def generate_altloc_selection(altloc_residue_dict, pdb):

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

    return altloc_selection

def from_altloc_generate_altloc_selection(pdb, altloc, occupancy_groups, params=master_phil.extract()):

    altloc_residue_dict = get_altloc_residue_dict(altloc, occupancy_groups)
    altloc_selection = generate_altloc_selection(altloc_residue_dict, pdb)

    return altloc_selection

def get_coincident_altloc_groups(altloc_groups_with_shared_residues, occupancy_groups, pdb):

    coincident_altloc_groups = []
    for altloc_group in altloc_groups_with_shared_residues:
        print("Altlocs with Shared residues {}".format(altloc_group))

        all_altloc_hier = []
        for altloc in altloc_group:
            altloc_residue_dict = get_altloc_residue_dict(altloc, occupancy_groups)
            altloc_hierarchy = generate_altloc_hiearchy(altloc_residue_dict, pdb)
            all_altloc_hier.append(altloc_hierarchy)

        # TODO Rewrite hierarchy overlap to handle any number of hierarchies

        if len(all_altloc_hier) > 2:
            raise Warning("Only two hierarchies can be compared currently. "
                          "Assuming that if a larger group, that they all overlap if the first two states overlap")

        if hierarchy_all_overlap(all_altloc_hier[0], all_altloc_hier[1], rmsd_distance=0.02):
            print("Altlocs {} are Coincident hierarchies".format(altloc_group))
            coincident_altloc_groups.append(altloc_group)

    return coincident_altloc_groups

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

# TODO Clean up how this sources the params
def get_coincident_altlocs(pdb, params=master_phil.extract()):

    occupancy_groups = get_occupancy_groups(pdb, params)
    altlocs = get_parameter_from_occupancy_groups(occupancy_groups, "altloc")
    altloc_groups_with_shared_residues = find_altlocs_with_shared_residues(altlocs, occupancy_groups)

    if len(altloc_groups_with_shared_residues) > 2:
        raise Warning("Current code can only handle two states")

    coincident_altloc_groups = get_coincident_altloc_groups(altloc_groups_with_shared_residues,
                                                            occupancy_groups, pdb)

    print ("Coincident_altloc_groups: {}".format(coincident_altloc_groups))

    return coincident_altloc_groups

# TODO Clean up how this sources the params
def get_altloc_hier(pdb, params=master_phil.extract()):
    coincident_altloc_groups = get_coincident_altlocs(pdb, params)
    occupancy_groups = get_occupancy_groups(pdb, params)

    for altloc_group in coincident_altloc_groups:
        altloc_residue_dict = get_altloc_residue_dict(altloc_group[0], occupancy_groups)
        altloc_hier = generate_altloc_hiearchy(altloc_residue_dict, pdb)

        yield altloc_hier

##############################################################
# From select_hierarchies.py
##############################################################

# from __future__ import print_function
# import sys
# import iotbx.pdb
# from scitbx.array_family import flex
# from giant.structure import calculate_paired_conformer_rmsds

def select_chains_to_vary(bound_pdb_path, ground_pdb_path):

    """ Select difference between two pdb files"""

    # Get protein hierarchies
    bound_pdb_inp = iotbx.pdb.input(bound_pdb_path)
    ground_pdb_inp = iotbx.pdb.input(ground_pdb_path)
    bound_ph = bound_pdb_inp.construct_hierarchy()
    ground_ph = ground_pdb_inp.construct_hierarchy()

    # TODO add option to work with other atoms (not just hetatms)
    bound_hetatm_hier = select_all_hetatm(bound_ph)
    ground_hetatm_hier = select_all_hetatm(ground_ph)
    bound_no_water_hier = remove_waters(bound_hetatm_hier)
    ground_no_water_hier = remove_waters(ground_hetatm_hier)

    # TODO Add a try excpet block and excpetions to catch when there are not differnces in both bound and ground
    bound_not_in_ground_chains, ground_not_in_bound_chains = compare_hier(bound_no_water_hier,ground_no_water_hier)

    bound_not_in_ground_hier, sel_bound_not_in_ground_hier = chains_select_hier(bound_not_in_ground_chains, bound_ph)
    ground_not_in_bound_hier, sel_ground_not_in_bound_hier = chains_select_hier(bound_not_in_ground_chains, ground_ph)

    # Are the differences between the bound and ground state within a rmsd cutoff
    if hierarchy_any_overlap(bound_not_in_ground_hier, ground_not_in_bound_hier, rmsd_distance=5):
        ground_bound_chains = list(set(bound_not_in_ground_chains + ground_not_in_bound_chains))
    else:
        ground_bound_chains = []

    # If there are common chains between bound not in ground and ground not in bound, check whether the common
    # chains lie within a rmsd cutoff of the chains found uniquely in the bound not in ground?

    #Is bound_ph appropriate here?

    common_chains = set(bound_not_in_ground_chains).intersection(ground_not_in_bound_chains)
    unique_bound_chains = set(bound_not_in_ground_chains).symmetric_difference(ground_not_in_bound_chains)

    common_chains_hier, common_chains_sel = chains_select_hier(common_chains, bound_ph)
    unique_bound_chains_hier, unique_bound_chains_sel = chains_select_hier(unique_bound_chains, bound_ph)

    if not hierarchy_any_overlap(common_chains_hier, unique_bound_chains_hier, rmsd_distance= 5):
        bound_not_in_ground_chains = list(set(bound_not_in_ground_chains).symmetric_difference(common_chains))
        ground_not_in_bound_chains = list(set(ground_not_in_bound_chains).symmetric_difference(common_chains))
        ground_bound_chains = list(set(bound_not_in_ground_chains + ground_not_in_bound_chains))

    return ground_bound_chains, bound_not_in_ground_chains, ground_not_in_bound_chains


def select_all_hetatm(protein_hier, cache=None, copy=True):
    if not cache: cache = protein_hier.atom_selection_cache()
    sel =cache.selection('hetero')
    return protein_hier.select(sel, copy_atoms = copy)


def remove_waters(protein_hier, cache=None, copy=True):
    if not cache: cache=protein_hier.atom_selection_cache()
    sel = cache.selection('(not resname HOH)')
    return protein_hier.select(sel, copy_atoms=copy)


def compare_hier(hier_1, hier_2, cache_1 = None, cache_2 = None, rmsd_cutoff = 0.1):

    """Take two protein hierarchies and find unique chains and chains with atoms that differ between the hierarchies"""

    # Generate caches for each hierarchy
    if not cache_1: cache_1 = hier_1.atom_selection_cache()
    if not cache_2: cache_2 = hier_2.atom_selection_cache()

    # If hierarchy is the same return None
    if hier_1.is_similar_hierarchy(hier_2):
        print ("Hierarchies are the same")
        return [], []

    # Get chains from hierarchies
    chain_ids_1 = [chain.id for chain in hier_1.only_model().chains()]
    chain_ids_2 = [chain.id for chain in hier_2.only_model().chains()]

    # Find common, and unique chains in each hierarchy
    common_chain_ids = list(set(chain_ids_1).intersection(chain_ids_2))

    chain_id1_not_in_id2 = list(set(chain_ids_1) - set(chain_ids_2))
    chain_id2_not_in_id1 = list(set(chain_ids_2) - set(chain_ids_1))

    # Find chains with atoms that are sperated by more than rmsd cutoff
    chains_with_distant_atoms = []
    for chain_id in common_chain_ids:
        hier_1_sel = cache_1.selection("chain {}".format(chain_id))
        hier_2_sel = cache_2.selection("chain {}".format(chain_id))

        hier_1_xyz = hier_1.atoms().extract_xyz().select(hier_1_sel)
        hier_2_xyz = hier_2.atoms().extract_xyz().select(hier_2_sel)

        if len(hier_1_xyz) == len(hier_2_xyz):

            distance = flex.mean((hier_1_xyz-hier_2_xyz).dot())

            if distance > rmsd_cutoff**2:
                print("Chain {} is different by {} A RMSD".format(chain_id, distance))
                chains_with_distant_atoms.append(chain_id)
        else:
            print ("Chain {} is a different length in ground and bound states.".format(chain_id) +
                   "Something has been removed?. Assuming chain is different")
            chains_with_distant_atoms.append(chain_id)

    # Add chains that have atoms in different positions to chains that are not shared between hierarchies
    chains_hier_1_different_heir_2 = chains_with_distant_atoms + chain_id1_not_in_id2
    chains_hier_2_different_heir_1 = chains_with_distant_atoms + chain_id2_not_in_id1

    return chains_hier_1_different_heir_2, chains_hier_2_different_heir_1


def chains_select_hier(chain_list, hier, cache = None, copy = True):

    """Return hierarchy and selection including only selected chains"""

    if not cache: cache=hier.atom_selection_cache()

    chains_sel = []
    for chain in chain_list:
        chains_sel.append("chain {}".format(chain))

    chains_selection_str = ' or '.join(chains_sel)
    sel = cache.selection(chains_selection_str)

    return hier.select(sel, copy_atoms=copy), sel


def hierarchy_any_overlap(hier_1, hier_2, cache_1 = None, cache_2 = None, rmsd_distance=5):

    """ Returns True if any atom in two hierarchies are within given rmsd distance"""

    all_distances = get_distance_between_hierarchies(hier_1, hier_2, cache_1, cache_2)

    overlap = False
    for distance in all_distances:
        if distance < rmsd_distance:
            overlap = True
            return overlap

    return overlap

def hierarchy_all_overlap(hier_1, hier_2, cache_1 = None, cache_2 = None, rmsd_distance=5):

    """ Returns True if all atom in two heirarchies are within given rmsd distance"""

    if not cache_1: cache_1 = hier_1.atom_selection_cache()
    if not cache_2: cache_2 = hier_2.atom_selection_cache()

    for chain_1, chain_2 in zip(hier_1.only_model().chains(),hier_2.only_model().chains()):
        print(type(chain_1), type(chain_2))
        for residue_group_1, residue_group_2 in zip(chain_1.residue_groups(), chain_2.residue_groups()):
            conformers_1 = residue_group_1.conformers()
            conformers_2 = residue_group_2.conformers()
            ret_list = calculate_paired_conformer_rmsds(conformers_1, conformers_2)
            rmsd_residue = ret_list[0][2]

            if rmsd_residue > rmsd_distance:
                print("RMSD {} > cutoff {} between chain {} resseq {} "
                      "and chain {} resseq {}".format(rmsd_residue, rmsd_distance,chain_1.id,residue_group_1.resseq,
                                                      chain_2.id,residue_group_1.resseq ))
                return False
    return True

def get_distance_between_hierarchies(hier_1, hier_2, cache_1 = None, cache_2 = None):

    """ Finds all distances between two hierarchies"""

    if not cache_1: cache_1 = hier_1.atom_selection_cache()
    if not cache_2: cache_2 = hier_2.atom_selection_cache()

    all_distances = []

    for atom_hier_1 in hier_1.atoms():
        for atom_hier_2 in hier_2.atoms():

            distance = atom_hier_1.distance(atom_hier_2)
            all_distances.append(distance)

    return all_distances
