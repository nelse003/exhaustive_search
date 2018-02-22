import os

from giant.structure.utils import transfer_residue_groups_from_other
from iotbx.pdb import hierarchy


def copy_atoms(path, prefix,start_xtal_num, end_xtal_num, new_ground_structure_path, atoms_new, out_dir):

    pdb_in = hierarchy.input(file_name=new_ground_structure_path)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()

    selection_string_list = []
    chains_new = set()
    for atom_new in atoms_new:
        selection_string = "(resid {} and chain {})".format(atom_new[1],atom_new[0])
        selection_string_list.append(selection_string)
        chains_new.add(atom_new[0])

    selection_string = "or".join(selection_string_list)
    new_atoms_sel = sel_cache.selection(selection_string)
    new_atoms_hier = pdb_in.hierarchy.select(new_atoms_sel)

    for num in range(start_xtal_num, end_xtal_num+1):
        xtal_name= prefix + "{0:0>4}".format(num)

        if os.path.exists(os.path.join(path,xtal_name,"refine.split.ground-state.pdb")):
            pdb_in_refine = hierarchy.input(file_name=os.path.join(path,xtal_name,"refine.split.ground-state.pdb"))
            acceptor_hierarchy = pdb_in_refine.construct_hierarchy()
            donor_hierarchy = new_atoms_hier
            acceptor_hier = transfer_residue_groups_from_other(acceptor_hierarchy, donor_hierarchy, in_place=False, verbose=False)

            if not os.path.exists(os.path.join(out_dir,xtal_name)):
                os.mkdir(os.path.join(out_dir,xtal_name))

            f = open(os.path.join(out_dir,xtal_name,"refine.split.ground-state_with_new_atoms.pdb"),"w")
            f.write(acceptor_hier.as_pdb_string(crystal_symmetry=pdb_in_refine.input.crystal_symmetry()))
            f.close()

            print(os.getcwd())

            os.chdir(os.path.join(out_dir,xtal_name))

            os.system('cp {} {}'.format(os.path.join(path,xtal_name,"refine.split.bound-state.pdb"),
                                        os.path.join(out_dir,xtal_name,"refine.split.bound-state.pdb")))


            os.system('giant.merge_conformations major={} minor={}'.format(
                os.path.join(out_dir,xtal_name,"refine.split.ground-state_with_new_atoms.pdb"),
                os.path.join(out_dir,xtal_name,"refine.split.bound-state.pdb")))

        else:
            print("pdb does not exist: {}".format(os.path.join(path,xtal_name,"refine.split.ground-state.pdb")))


copy_atoms(path = "/dls/labxchem/data/2018/lb18145-55/processing/analysis/initial_model",
           prefix = "NUDT22A-x",
           start_xtal_num=421,
           end_xtal_num=421,
           new_ground_structure_path="/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/"\
                                     "NUDT22A_new_ground_state_x1058/ground_state_refine_2_more_waters.pdb",
           atoms_new=[['W','230'],['W','231'],['W','232'],['W','233'],['W','234'],['W','235']],
           out_dir="/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/occupancy_group_test")