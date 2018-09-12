import os

from giant.structure.utils import transfer_residue_groups_from_other
from iotbx.pdb import hierarchy


def copy_atoms(path, prefix, start_xtal_num, end_xtal_num,
               new_ground_structure_path, atoms_new, out_dir):

    """ Copy one set of atoms a pdb structure to other pdb structures
    Also run refinement
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

    # xtals = ['NUDT22A-x0243', 'NUDT22A-x0421','NUDT22A-x0391']
    xtals = []
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    for xtal_name in xtals:

        if os.path.exists(os.path.join(path, xtal_name,
                                       "refine.split.ground-state.pdb")):

            pdb_in_refine = hierarchy.input(
                file_name=os.path.join(path, xtal_name,
                                       "refine.split.ground-state.pdb"))

            acceptor_hierarchy = pdb_in_refine.construct_hierarchy()
            donor_hierarchy = new_atoms_hier
            acceptor_hier = transfer_residue_groups_from_other(
                acceptor_hierarchy, donor_hierarchy, in_place=False,
                verbose=False)

            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            if not os.path.exists(os.path.join(out_dir, xtal_name)):
                os.mkdir(os.path.join(out_dir, xtal_name))

            f = open(
                os.path.join(out_dir,
                             xtal_name,
                             "refine.split.ground-state_with_new_atoms.pdb"),
                "w")

            f.write(acceptor_hier.as_pdb_string(
                crystal_symmetry=pdb_in_refine.input.crystal_symmetry()))
            f.close()

            print(os.getcwd())
            os.chdir(os.path.join(out_dir, xtal_name))
            os.system('cp {} {}'.format(
                os.path.join(path, xtal_name, "refine.split.bound-state.pdb"),
                os.path.join(out_dir, xtal_name,
                             "refine.split.bound-state.pdb")))

            os.system('cp {} {}'.format(os.path.join(path, xtal_name, "*.cif"),
                                        os.path.join(out_dir, xtal_name)))

            os.system('cp -rL {} {}'.format(os.path.join(path, xtal_name,
                                                         "refine.mtz"),
                                        os.path.join(out_dir, xtal_name,
                                                     "refine.mtz")))

            os.system('giant.merge_conformations major={} minor={}'.format(
                os.path.join(out_dir, xtal_name,
                             "refine.split.ground-state_with_new_atoms.pdb"),
                os.path.join(out_dir, xtal_name,
                             "refine.split.bound-state.pdb")))

            cmds = "source /dls/science/groups/i04-1/software/" \
            "pandda-update/ccp4/ccp4-7.0/setup-scripts/ccp4.setup-sh \n"

            cmds += "giant.quick_refine {} {} {} params={} \n".format(
                os.path.join(out_dir, xtal_name, "multi-state-model.pdb"),
                os.path.join(out_dir, xtal_name, "refine.mtz"),
                os.path.join(out_dir, xtal_name, "*.cif"),
                os.path.join(out_dir, xtal_name,
                             "multi-state-restraints.refmac.params"))

            os.system(cmds)

        else:
            print("pdb does not exist: {}".format(
                os.path.join(path, xtal_name,
                             "refine.split.ground-state.pdb")))


def copy_covalent_ratios(path, prefix, start_xtal_num, end_xtal_num,
                         new_ground_structure_path, atoms_new, atoms_remove, out_dir, qsub):

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

        # for quick rerun
        # if os.path.exists(os.path.join(out_dir,xtal_name,"refine.pdb")):
        #     continue

        if os.path.exists(os.path.join(path, xtal_name,
                                       "dimple.pdb")):

            pdb_in_refine = hierarchy.input(
                file_name=os.path.join(path, xtal_name,
                                       "dimple.pdb"))


            acceptor_hierarchy = pdb_in_refine.hierarchy()
            #remove atoms
            refine_sel_cache = pdb_in_refine.hierarchy.atom_selection_cache()

            print(not_selection_string)

            remove_atoms_sel = refine_sel_cache.selection(not_selection_string)
            removed_hier = acceptor_hierarchy.select(remove_atoms_sel)

            # Add atoms
            donor_hierarchy = new_atoms_hier
            acceptor_hier = transfer_residue_groups_from_other(
                removed_hier, donor_hierarchy, in_place=False,
                verbose=False)


            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            if not os.path.exists(os.path.join(out_dir, xtal_name)):
                os.mkdir(os.path.join(out_dir, xtal_name))

            f = open(
                os.path.join(out_dir,
                             xtal_name,
                             "dimple_with_lig.pdb"),
                "w+")

            #f.write("LINKR        C   LIG E   1                 SG ACYS A  73                LIG-CYS\n")
            f.write(acceptor_hier.as_pdb_string(
                crystal_symmetry=pdb_in_refine.input.crystal_symmetry()))


            # # Add Link record to dimple_with_lig
            # for line in f:
            #     if line.startswith("HEADER"):
            #         continue
            #     if line.startswith("COMPND"):
            #         continue
            #     if line.startswith("REMARK"):
            #         continue
            #     else:
            #         f.write("LINKR        C   LIG E   1                 SG ACYS A  73                LIG-CYS")
            #         break

            f.close()


            print(os.getcwd())
            os.chdir(os.path.join(out_dir, xtal_name))
            os.system('cp {} {}'.format(
                os.path.join(path, xtal_name, "dimple.pdb"),
                os.path.join(out_dir, xtal_name,
                             "dimple.pdb")))

            os.system('cp {} {}'.format(os.path.join("/dls/science/groups/i04-1/elliot-dev/Work/"
                                                     "exhaustive_search_data/NUDT7_covalent/NUDT7A-x1812/"
                                                     "NUDT7A-x1812LIG-CYS.cif"),
                                        os.path.join(out_dir, xtal_name, "{}_LIG_CYS.cif".format(xtal_name))))

            print(os.listdir(os.path.join(out_dir, xtal_name)))

            os.system('cp -rL {} {}'.format(os.path.join(path, xtal_name,
                                                         "dimple.mtz"),
                                        os.path.join(out_dir, xtal_name,
                                                     "dimple.mtz")))

            os.system('giant.merge_conformations major={} minor={}'.format(
                os.path.join(out_dir, xtal_name,
                             "dimple.pdb"),
                os.path.join(out_dir, xtal_name,
                             "dimple_with_lig.pdb")))

            with open(os.path.join(out_dir, xtal_name,"multi-state-model.pdb"),"r") as original:
                multi_model = original.read()
            with open(os.path.join(out_dir, xtal_name,"multi-state-model.pdb"),"w") as modified:
                modified.write("LINKR        C  CLIG E   1                 SG ACYS A  73                LIG-CYS\n")
                modified.write("LINKR        C  DLIG E   1                 SG ACYS A  73                LIG-CYS\n")
                modified.write(multi_model)

            cmds = "source /dls/science/groups/i04-1/software/" \
            "pandda-update/ccp4/ccp4-7.0/setup-scripts/ccp4.setup-sh \n"

            cmds += "giant.quick_refine {} {} {} params={} \n".format(
                os.path.join(out_dir, xtal_name, "multi-state-model.pdb"),
                os.path.join(out_dir, xtal_name, "dimple.mtz"),
                os.path.join(out_dir, xtal_name, "*.cif"),
                os.path.join(out_dir, xtal_name,
                             "multi-state-restraints.refmac.params"))
            if qsub:
                f = open(
                    os.path.join(out_dir,
                                 xtal_name,
                                 "{}_quick_refine.sh".format(xtal_name)),
                    "w")

                f.write(cmds)
                f.close()

                os.system('qsub {}'.format(os.path.join(out_dir,xtal_name,"{}_quick_refine.sh".format(xtal_name))))
            else:
                os.system(cmds)

        else:
            print("pdb does not exist: {}".format(
                os.path.join(path, xtal_name,
                             "dimple.pdb")))


def copy_titration(path, prefix, start_xtal_num, end_xtal_num,
                         new_ground_structure_path, atoms_new, atoms_remove, out_dir, qsub):
    """ Script that works only with titration data for initial refinement 

    Copy dimple pdb, mtz and cif 
    Copy ligand atoms from existing coordinates 
    Remove water atoms
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

        # for quick rerun
        if os.path.exists(os.path.join(out_dir, xtal_name, "refine.pdb")):
            continue

        if os.path.exists(os.path.join(path, xtal_name,
                                       "dimple.pdb")):

            pdb_in_refine = hierarchy.input(
                file_name=os.path.join(path, xtal_name,
                                       "dimple.pdb"))

            acceptor_hierarchy = pdb_in_refine.construct_hierarchy()

            #remove atoms
            refine_sel_cache = pdb_in_refine.hierarchy.atom_selection_cache()

            print(not_selection_string)

            remove_atoms_sel = refine_sel_cache.selection(not_selection_string)
            removed_hier = acceptor_hierarchy.select(remove_atoms_sel)

            # add atoms
            donor_hierarchy = new_atoms_hier
            acceptor_hier = transfer_residue_groups_from_other(
                removed_hier, donor_hierarchy, in_place=False,
                verbose=False)


            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

            if not os.path.exists(os.path.join(out_dir, xtal_name)):
                os.mkdir(os.path.join(out_dir, xtal_name))

            f = open(
                os.path.join(out_dir,
                             xtal_name,
                             "dimple_with_lig.pdb"),
                "w+")

            f.write(acceptor_hier.as_pdb_string(
                crystal_symmetry=pdb_in_refine.input.crystal_symmetry()))

            f.close()

            print(os.getcwd())
            os.chdir(os.path.join(out_dir, xtal_name))
            os.system('cp {} {}'.format(
                os.path.join(path, xtal_name, "dimple.pdb"),
                os.path.join(out_dir, xtal_name,
                             "dimple.pdb")))

            # os.system('cp {} {}'.format(os.path.join("/dls/labxchem/data/2017/lb18145-3/processing/"
            #                                          "analysis/initial_model/NUDT7A-x1237/NUOOA000181a.cif"),
            #                             os.path.join(out_dir, xtal_name, "NUOOA000181a.cif".format(xtal_name))))

            os.system('cp {} {}'.format(os.path.join("/dls/science/groups/i04-1/elliot-dev/Work/"
                                                     "exhaustive_search_data/NUDT7_covalent/NUDT7A-x1812/"
                                                     "NUDT7A-x1812LIG-CYS.cif"),
                                        os.path.join(out_dir, xtal_name, "{}_LIG_CYS.cif".format(xtal_name))))

            print(os.listdir(os.path.join(out_dir, xtal_name)))

            os.system('cp -rL {} {}'.format(os.path.join(path, xtal_name,
                                                         "dimple.mtz"),
                                            os.path.join(out_dir, xtal_name,
                                                         "dimple.mtz")))

            os.system('giant.merge_conformations major={} minor={}'.format(
                os.path.join(out_dir, xtal_name,
                             "dimple.pdb"),
                os.path.join(out_dir, xtal_name,
                             "dimple_with_lig.pdb")))


            cmds = "source /dls/science/groups/i04-1/software/" \
                   "pandda-update/ccp4/ccp4-7.0/setup-scripts/ccp4.setup-sh \n"

            cmds += "giant.quick_refine {} {} {} params={} \n".format(
                os.path.join(out_dir, xtal_name, "multi-state-model.pdb"),
                os.path.join(out_dir, xtal_name, "dimple.mtz"),
                os.path.join(out_dir, xtal_name, "*.cif"),
                os.path.join(out_dir, xtal_name,
                             "multi-state-restraints.refmac.params"))
            if qsub:
                f = open(
                    os.path.join(out_dir,
                                 xtal_name,
                                 "{}_quick_refine.sh".format(xtal_name)),
                    "w")

                f.write(cmds)
                f.close()

                os.system('qsub {}'.format(os.path.join(out_dir, xtal_name, "{}_quick_refine.sh".format(xtal_name))))
            else:
                os.system(cmds)

        else:
            print("pdb does not exist: {}".format(
                os.path.join(path, xtal_name,
                             "dimple.pdb")))


# copy_titration(path="/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model",
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

copy_covalent_ratios(path="/dls/labxchem/data/2017/lb18145-68/processing/initial_model",
                     prefix='NUDT7A-x',
                     start_xtal_num=6192,
                     end_xtal_num=6251,
                     new_ground_structure_path="/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/NUDT7_covalent/NUDT7A-x1812/refine.pdb",
                     atoms_new=[['E','1'],
                     atoms_remove = [['A','196']]
                     out_dir="/dls/science/groups/i04-1/elliot-dev/Work/"
                   "exhaustive_search_data/covalent_ratios_dose",
                     qsub = True)

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