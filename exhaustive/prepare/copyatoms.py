import os

from giant.structure.utils import transfer_residue_groups_from_other
from iotbx.pdb import hierarchy
from copy_phil import copy_phil

def copy_atoms(copy_params):

    """ Copy atoms from one pdb file to many, then refine.
    
    Copy dimple pdb, mtz and cif with cys bond
    Copy ligand atoms from existing coordinates 
    Run giant.merge_conformations to generate a multi state model
    Copy link records suitable for both conformers of the ligand
    Run quick refine to generate refined ligand 
    """,

    # generate output directory if it doesn't exist
    if not os.path.exists(copy_params.output.out_dir):
        os.mkdir(copy_params.output.out_dir)

    # read in PDB file from which atoms are to be taken from (ground structure)
    pdb_in = hierarchy.input(file_name=copy_params.input.base_pdb)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()

    # produce a hierarchy with atoms to copied
    selection_string_list = []
    chains_new = set()
    for atom_new in copy_params.input.atoms_new:
        selection_string = "(resid {} and chain {})".format(atom_new[1],
                                                            atom_new[0])
        selection_string_list.append(selection_string)
        chains_new.add(atom_new[0])
    selection_string = "or".join(selection_string_list)
    new_atoms_sel = sel_cache.selection(selection_string)
    new_atoms_hier = pdb_in.hierarchy.select(new_atoms_sel)

    # Produce a selection string to determine which atoms are removed
    selection_string_list = []
    if copy_params.input.atoms_remove is not None:
        for atom_remove in copy_params.input.atoms_remove:
            selection_string = "(resid {} and chain {})".format(atom_remove[1],
                                                                atom_remove[0])
            selection_string_list.append(selection_string)

        selection_string = "or".join(selection_string_list)
        not_selection_string ="not ({})".format(selection_string)

    # Define xtals to loop over
    xtals = copy_params.input.xtal_list
    for num in range(copy_params.input.start_xtal_number,
                     copy_params.input.end_xtal_number + 1):

        xtal_name = copy_params.input.prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    # Loop over all xtals
    for xtal_name in xtals:

        # For quick rerun
        if os.path.exists(os.path.join(copy_params.output.out_dir,
                                       xtal_name,
                                       copy_params.output.refine_pdb)) and \
            not copy_params.settings.overwrite:
            continue

        # Run only if sufficent input data
        if not os.path.exists(os.path.join(copy_params.input.path,
                                           xtal_name,
                                           copy_params.input.pdb_style)):

            print("pdb does not exist: {}".format(
                os.path.join(copy_params.input.path,
                             xtal_name, copy_params.input.pdb_style)))
            continue

        pdb_in_refine = hierarchy.input(
            file_name=os.path.join(copy_params.input.path, xtal_name, copy_params.input.pdb_style))

        acceptor_hierarchy = pdb_in_refine.construct_hierarchy()

        #remove atoms from xtal
        if copy_params.input.atoms_remove is not None:
            refine_sel_cache = pdb_in_refine.hierarchy.atom_selection_cache()
            remove_atoms_sel = refine_sel_cache.selection(not_selection_string)
            removed_hier = acceptor_hierarchy.select(remove_atoms_sel)
            working_hier = removed_hier
        else:
            working_hier = acceptor_hierarchy

        # Add atoms from base_pdb
        donor_hierarchy = new_atoms_hier
        acceptor_hier = transfer_residue_groups_from_other(
            working_hier, donor_hierarchy, in_place=False,
            verbose=False)

        # Generate output xtal directories
        if not os.path.exists(os.path.join(
                copy_params.output.out_dir, xtal_name)):
            os.mkdir(os.path.join(copy_params.output.out_dir, xtal_name))

        # Write output pdb with changed atoms
        f = open(os.path.join(copy_params.output.out_dir,xtal_name,
                              copy_params.output.pdb),"w+")
        f.write(acceptor_hier.as_pdb_string(
            crystal_symmetry=pdb_in_refine.input.crystal_symmetry()))
        f.close()

        # Copy the input pdb to output directory
        os.chdir(os.path.join(copy_params.output.out_dir, xtal_name))
        os.system('cp {} {}'.format(
            os.path.join(copy_params.input.path,
                         xtal_name, copy_params.input.pdb_style),
            os.path.join(copy_params.output.out_dir,
                         xtal_name, copy_params.input.pdb_style)))

        # Copy the input cif to output_directory
        os.system('cp {} {}'.format(copy_params.input.cif,
                                    os.path.join(copy_params.output.out_dir,
                                                 xtal_name,
                                                 os.path.basename(copy_params.input.cif))))

        # Copy the input mtz to output directory
        os.system('cp -rL {} {}'.format(os.path.join(copy_params.input.pdb,
                                                     xtal_name,
                                                     copy_params.input.mtz_style),
                                    os.path.join(copy_params.output.out_dir,
                                                 xtal_name,
                                                 copy_params.input.mtz_style)))
        # Run giant.merge_conforamtions
        os.system('giant.merge_conformations major={} minor={}'.format(
            os.path.join(copy_params.output.out_dir,
                         xtal_name, copy_params.input.pdb_style),
            os.path.join(copy_params.output.out_dir,
                         xtal_name, copy_params.output.pdb)))


        # Add link record strings into multimodel pdb file, prior to refinement
        if copy_params.input.link_record_list is not None:

            with open(os.path.join(
                    copy_params.output.out_dir,
                    xtal_name,
                    copy_params.output.multi_state_model_pdb),
                    "r") as original:

                multi_model = original.read()

            with open(os.path.join(
                    copy_params.output.out_dir,
                    xtal_name,
                    copy_params.output.multi_state_model_pdb),"w") as modified:

                for link_record in copy_params.input.link_record_list:
                    modified.write(link_record)

                modified.write(multi_model)

        # Run giant.quick_refine
        cmds = "source {}".format(copy_params.settings.ccp4_path)

        cmds += "giant.quick_refine {} {} {} params={} \n".format(
            os.path.join(copy_params.output.out_dir,
                         xtal_name,
                         copy_params.output.multi_state_model_pdb),
            os.path.join(copy_params.output.out_dir,
                         xtal_name, copy_params.input.mtz_style),
            os.path.join(copy_params.output.out_dir,
                         xtal_name,
                         copy_params.input.cif),
            os.path.join(copy_params.output.out_dir,
                         xtal_name, copy_params.settings.param))

        if copy_params.settings.qsub:
            f = open(
                os.path.join(copy_params.output.out_dir,xtal_name,
                             "{}_quick_refine.sh".format(xtal_name)),"w")

            f.write(cmds)
            f.close()

            os.system('qsub {}'.format(os.path.join(copy_params.output.out_dir,
                                                    xtal_name,
                                                    "{}_quick_refine.sh".format(xtal_name))))
        else:
            os.system(cmds)


#titration
# (path="/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model",
#                      prefix='NUDT7A-x',
#                      start_xtal_num=2006,
#                      end_xtal_num=2073,
#                      new_ground_structure_path="/dls/labxchem/data/2017/lb18145-3/processing/analysis/initial_model/NUDT7A-x1237/refine.pdb",
#                      atoms_new=[['E','1']],
#                      atoms_remove=[['B','60'],['B','151'],['B','189'],['B','33'],['B','11'],['B','40'],['B','196']],
#                      out_dir="/dls/science/groups/i04-1/elliot-dev/Work/"
#                    "exhaustive_search_data/titration_series",
#                      qsub = True)

# To be used for covalent atoms

copy_params = copy_phil.extract()

copy_params.input.path = "/dls/labxchem/data/2018/lb18145-68/processing/analysis/initial_model"
copy_params.input.prefix = 'NUDT7A-x'
copy_params.input.start_xtal_number = 6192
copy_params.input.end_xtal_number = 6251
copy_params.input.base_pdb = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/" \
                             "NUDT7_covalent/NUDT7A-x1812/refine.pdb"
copy_params.input.atoms_new = [['E','1']]
copy_params.input.atoms_remove = [['B','196']]
copy_params.input.cif = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                        "exhaustive_search_data/NUDT7_covalent" \
                        "/NUDT7A-x1812/NUDT7A-x1812LIG-CYS.cif"
copy_params.input.link_record_list =["LINKR        C  CLIG E   1                 SG ACYS A  73                LIG-CYS\n",
                                     "LINKR        D  CLIG E   1                 SG ACYS A  73                LIG-CYS\n"]

copy_params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                   "exhaustive_search_data/test_copy_atoms"

copy_params.settings.overwrite = True

copy_atoms(copy_params)

# Commented out for testing of new dimple based function
# copy_atoms(path="/dls/labxchem/data/2018/lb18145-55/processing/analysis/"
#                 "initial_model",
#            prefix="NUDT22A-x",
#            start_xtal_num=182,
#            end_xtal_num=182,
#            new_ground_structure_path="/dls/science/groups/i04-1/elliot-dev/"
#                                      "Work/exhaustive_search_data/"
#                                      "NUDT22A_new_ground_state_x1058/"
#                                      "ground_state_refine_2_more_waters.pdb",
#            atoms_new=[['W', '230'], ['W', '231'], ['W', '232'], ['W', '233'],
#                       ['W', '234'], ['W', '235']],
#            out_dir="/dls/science/groups/i04-1/elliot-dev/Work/"
#                    "exhaustive_search_data/occupancy_group_with_refinement")

# NUDT22 [['W','11'],['W','230'],['W','231'],['W','232'],
# ['W','233'],['W','234'],['W','235']]