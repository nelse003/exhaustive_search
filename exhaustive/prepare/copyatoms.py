import os

from giant.structure.utils import transfer_residue_groups_from_other
from iotbx.pdb import hierarchy

def copy_covalent_ratios(path, prefix, start_xtal_num, end_xtal_num,
                         new_ground_structure_path, atoms_new,
                         atoms_remove, out_dir, qsub,
                         overwrite, input_pdb, output_pdb, refine_pdb,
                         input_cif, output_cif_suffix, input_mtz,
                         multi_state_model_pdb, link_record_list,
                         ccp4_path, param_file):

    """ Script that works only with covalent ratio data for initial refinement 
    
    Copy dimple pdb, mtz and cif with cys bond
    Copy ligand atoms from existing coordinates (NUDT7A-x1812)
    Run giant.merge_conformations to generate a multi state model
    Copy link records suitable for both conformers of the ligand
    Run quick refine to generate refined ligand 
    
    Needs to be trimmed and turned into reusable code
    """,

    pdb_in = hierarchy.input(file_name=new_ground_structure_path)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()

    selection_string_list = []
    chains_new = set()
    for atom_new in atoms_new:
        selection_string = "(resid {} and chain {})".format(atom_new[1],
                                                            atom_new[0])
        selection_string_list.append(selection_string)
        chains_new.add(atom_new[0])

    selection_string = "or".join(selection_string_list)
    new_atoms_sel = sel_cache.selection(selection_string)
    new_atoms_hier = pdb_in.hierarchy.select(new_atoms_sel)

    selection_string_list = []

    if atoms_remove is not None:
        for atom_remove in atoms_remove:
            selection_string = "(resid {} and chain {})".format(atom_remove[1],
                                                                atom_remove[0])
            selection_string_list.append(selection_string)

        selection_string = "or".join(selection_string_list)
        not_selection_string ="not ({})".format(selection_string)

    xtals = []
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    for xtal_name in xtals:

        #for quick rerun
        if os.path.exists(os.path.join(out_dir, xtal_name, refine_pdb)) and \
            not overwrite:
            continue

        if not os.path.exists(os.path.join(path, xtal_name,input_pdb)):
            print("pdb does not exist: {}".format(os.path.join(path, xtal_name, input_pdb)))

        else:

            pdb_in_refine = hierarchy.input(
                file_name=os.path.join(path, xtal_name, input_pdb))

            acceptor_hierarchy = pdb_in_refine.construct_hierarchy()

            #remove atoms
            if atoms_remove is not None:
                refine_sel_cache = pdb_in_refine.hierarchy.atom_selection_cache()
                remove_atoms_sel = refine_sel_cache.selection(not_selection_string)
                removed_hier = acceptor_hierarchy.select(remove_atoms_sel)
                working_hier = removed_hier
            else:
                working_hier = acceptor_hierarchy

            # Add atoms
            donor_hierarchy = new_atoms_hier
            acceptor_hier = transfer_residue_groups_from_other(
                working_hier, donor_hierarchy, in_place=False,
                verbose=False)

            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            if not os.path.exists(os.path.join(out_dir, xtal_name)):
                os.mkdir(os.path.join(out_dir, xtal_name))

            f = open(os.path.join(out_dir,xtal_name, output_pdb),"w+")

            f.write(acceptor_hier.as_pdb_string(
                crystal_symmetry=pdb_in_refine.input.crystal_symmetry()))
            f.close()

            os.chdir(os.path.join(out_dir, xtal_name))
            os.system('cp {} {}'.format(
                os.path.join(path, xtal_name, input_pdb),
                os.path.join(out_dir, xtal_name, input_pdb)))

            os.system('cp {} {}'.format(input_cif, os.path.join(out_dir,
                    xtal_name, "{}_{}".format(xtal_name,
                                              output_cif_suffix))))

            os.system('cp -rL {} {}'.format(os.path.join(path, xtal_name,
                                                         input_mtz),
                                        os.path.join(out_dir, xtal_name,
                                                     input_mtz)))

            os.system('giant.merge_conformations major={} minor={}'.format(
                os.path.join(out_dir, xtal_name, input_pdb),
                os.path.join(out_dir, xtal_name, output_pdb)))

            if link_record_list is not None:
                with open(os.path.join(out_dir, xtal_name, multi_state_model_pdb),"r") as original:
                    multi_model = original.read()
                with open(os.path.join(out_dir, xtal_name, multi_state_model_pdb),"w") as modified:
                    for link_record in link_record_list:
                        modified.write(link_record)
                    modified.write(multi_model)

            cmds = "source {}".format(ccp4_path)

            cmds += "giant.quick_refine {} {} {} params={} \n".format(
                os.path.join(out_dir, xtal_name, multi_state_model_pdb),
                os.path.join(out_dir, xtal_name, input_mtz),
                os.path.join(out_dir, xtal_name, "{}_{}".format(xtal_name, output_cif_suffix)),
                os.path.join(out_dir, xtal_name, param_file))
            if qsub:
                f = open(
                    os.path.join(out_dir,xtal_name,
                                 "{}_quick_refine.sh".format(xtal_name)),"w")

                f.write(cmds)
                f.close()

                os.system('qsub {}'.format(os.path.join(out_dir,xtal_name,"{}_quick_refine.sh".format(xtal_name))))
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

copy_covalent_ratios(path="/dls/labxchem/data/2018/lb18145-68/processing/analysis/initial_model",
                     prefix='NUDT7A-x',
                     start_xtal_num=6192,
                     end_xtal_num=6251,
                     new_ground_structure_path="/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/NUDT7_covalent/NUDT7A-x1812/refine.pdb",
                     atoms_new=[['E','1']],
                     atoms_remove = [['B','196']],
                     out_dir="/dls/science/groups/i04-1/elliot-dev/Work/"
                   "exhaustive_search_data/test_copy_atoms",
                     qsub = True,
                     overwrite = True,
                     input_pdb = "dimple.pdb",
                     output_pdb = "dimple_with_lig.pdb",
                     refine_pdb = "refine.pdb",
                     input_cif = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                                 "exhaustive_search_data/NUDT7_covalent" \
                                 "/NUDT7A-x1812/NUDT7A-x1812LIG-CYS.cif",
                     output_cif_suffix = "LIG_CYS.cif",
                     input_mtz = "dimple.mtz",
                     multi_state_model_pdb = "multi-state-model.pdb",
                     link_record_list = ["LINKR        C  CLIG E   1                 SG ACYS A  73                LIG-CYS\n",
                                         "LINKR        D  CLIG E   1                 SG ACYS A  73                LIG-CYS\n"],
                     ccp4_path="/dls/science/groups/i04-1/software/" \
            "pandda-update/ccp4/ccp4-7.0/setup-scripts/ccp4.setup-sh \n",
                     param_file = "multi-state-restraints.refmac.params")

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