import os
import sys
import shutil
from collections import Counter

sys.path.append("/dls/science/groups/i04-1/elliot-dev/parse_xchemdb")
from refinement.merge_cif import merge_cif
from refinement.prepare_scripts import write_refmac_csh
from refinement.prepare_scripts import write_exhaustive_csh
from refinement.prepare_scripts import write_phenix_csh
from refinement.prepare_scripts import write_buster_csh

if __name__ == "__main__":
    """ Submit qsub refinements for all repeat soaks
    
    Requires source envs/parse_xchemdb
    """
    xtals_with_compound = {}

    # DCP2B FMOPL00213a
    # What are these hits?

    # # # DCP2B FMOPL000435a
    # xtals_with_compound["DCP2B-x0146"] = 'FMOPL000435a'
    # start_xtal_num = 993
    # end_xtal_num = 1048
    # for num in range(start_xtal_num, end_xtal_num + 1):
    #     xtal_name = "DCP2B-x" + "{0:0>4}".format(num)
    #     xtals_with_compound[xtal_name] = 'FMOPL000435a'
    #
    # # NUDT7A OX210 1774-1810, 0299
    # xtals_with_compound['NUDT7A-x0299'] = 'OX210'
    # start_xtal_num = 1774
    # end_xtal_num = 1810
    # for num in range(start_xtal_num, end_xtal_num + 1):
    #     xtal_name = "NUDT7A-x" + "{0:0>4}".format(num)
    #     xtals_with_compound[xtal_name] = 'OX210'
    #
    # # NUDT7A NUOOA000181a
    # # 1739 - 1773
    # # 2006 - 2073
    # start_xtal_num = 1739
    # end_xtal_num = 1773
    # for num in range(start_xtal_num, end_xtal_num + 1):
    #     xtal_name = "NUDT7A-x" + "{0:0>4}".format(num)
    #     xtals_with_compound[xtal_name] = 'NUOOA000181a'
    # start_xtal_num = 2006
    # end_xtal_num = 2073
    # for num in range(start_xtal_num, end_xtal_num + 1):
    #     xtal_name = "NUDT7A-x" + "{0:0>4}".format(num)
    #     xtals_with_compound[xtal_name] = 'NUOOA000181a'
    #
    #
    # # # NUDT22A FMOPL00622a
    # # from DSI poised library (x0938-x0976)
    # xtals_with_compound['NUDT22A-x0182'] = 'FMOPL00622a_DSPL'
    # start_xtal_num = 909
    # end_xtal_num = 937
    # for num in range(start_xtal_num, end_xtal_num + 1):
    #     xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
    #     xtals_with_compound[xtal_name] = 'FMOPL00622a_DSPL'
    #
    # # # NUDT22A FMOPL00622a
    # # from DSI poised library (x0938-x0976)
    # start_xtal_num = 938
    # end_xtal_num = 976
    # for num in range(start_xtal_num, end_xtal_num + 1):
    #     xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
    #     xtals_with_compound[xtal_name] = 'FMOPL000622a_DSI_poised'
    #
    # # NUDT22A 133725a x0421, x1040 - x1059
    # start_xtal_num = 1040
    # end_xtal_num = 1059
    # xtals_with_compound["NUDT22A-x0421"] = "N13725a"
    # for num in range(start_xtal_num, end_xtal_num + 1):
    #     xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
    #     xtals_with_compound[xtal_name] = "N13725a"
    #
    # # NUDT22A 13369a x0243, x977 to x1008
    #
    # start_xtal_num = 977
    # end_xtal_num = 1008
    # xtals_with_compound['NUDT22A-x0243'] = 'N13369a'
    # for num in range(start_xtal_num, end_xtal_num + 1):
    #     xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
    #     xtals_with_compound[xtal_name] = 'N13369a'

    #NUDT22A 13663a x0391, x1009 to x1039
    xtals_with_compound['NUDT22A-x0391'] = 'N13663a'
    start_xtal_num = 1009
    end_xtal_num = 1039
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = "N13663a"

    compound_codes = set()
    for compounds in xtals_with_compound.values():
        compound_codes.add(compounds)

    for xtal, compound in xtals_with_compound.items():

        if "DCP2B" in xtal:
            data_folder = (
                "/dls/science/groups/i04-1/elliot-dev/Work/"
                "exhaustive_search_data/DCP2B_18_09_20_exhaus"
            )
            input_cif = (
                "/dls/labxchem/data/2016/lb13385-64/processing/"
                "analysis/initial_model/DCP2B-x0146/FMOPL000435a.cif"
            )

        elif "NUDT22" in xtal:
            print(compound)
            data_folder_pre = (
                "/dls/science/groups/i04-1/elliot-dev/Work/"
                "exhaustive_search_data/NUDT22_repeat_soaks"
            )

            if compound == "FMOPL000622a_DSI_poised":
                folder_ext = "FMOPL00622a_Alternate_Diastereomer"
                input_cif = (
                    "/dls/labxchem/data/2018/lb18145-55/processing/"
                    "analysis/initial_model/NUDT22A-x0938/"
                    "FMOPL00622a_AlternateDiastereomer.cif"
                )

            elif compound == "FMOPL00622a_DSPL":
                folder_ext = "FMOPL00622a"
                input_cif = (
                    "/dls/labxchem/data/2018/lb18145-55/processing/"
                    "analysis/initial_model/NUDT22A-x0937/"
                    "FMOPL00622a.cif"
                )

            elif compound == "N13725a":
                folder_ext = compound
                input_cif = (
                    "/dls/labxchem/data/2018/lb18145-55/"
                    "processing/analysis/initial_model/NUDT22A-x0421/N13725a.cif"
                )

            elif compound == "N13369a":
                folder_ext = compound
                input_cif = (
                    "/dls/labxchem/data/2018/lb18145-55/"
                    "processing/analysis/initial_model/NUDT22A-x0977/"
                    "N13369a.cif"
                )

            elif compound == "N13663a":
                folder_ext = compound
                input_cif = (
                    "/dls/labxchem/data/2018/lb18145-55/"
                    "processing/analysis/initial_model/NUDT22A-x0391/"
                    "N13663a.cif"
                )

            data_folder = os.path.join(data_folder_pre, folder_ext)

        elif "NUDT7" in xtal:
            if compound == "OX210":
                data_folder = (
                    "/dls/science/groups/i04-1/elliot-dev/Work/"
                    "exhaustive_search_data/NUDT7_Copied_atoms/OX210"
                )

                input_cif = (
                    "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"
                    "NUDT7_Copied_atoms/OX210/NUDT7A-x1787/OX-210.cif"
                )
            else:
                # This isn't copied atoms. Can't find that...
                data_folder = "/dls/labxchem/data/2017/lb18145-3/processing/analysis/initial_model"
                input_cif = (
                    "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model"
                    "/NUDT7A-x1739/NUOOA000181a.cif"
                )

        # cif file issue management
        issue_compounds =["FMOPL00622a_DSPL",
                           "FMOPL000622a_DSI_poised",
                          "FMOPL000435a",
                          "N13369a",
                          "N13663a"]

        input_pdb = os.path.join(data_folder, xtal, "refine.pdb")
        input_mtz = os.path.join(data_folder, xtal, "refine.mtz")

        if not os.path.isdir(os.path.join(data_folder, xtal)):
            continue

        print("HERE")

        if not os.path.exists(input_pdb):
            for pdb_file in os.listdir(os.path.join(data_folder, xtal)):
                if "refine" in pdb_file and ".pdb" in pdb_file:
                    input_pdb = os.path.join(data_folder, xtal, pdb_file)
                else:
                    print("PDB doesn't exist, skipping: {}".format(input_pdb))
                    continue

        if not os.path.exists(input_mtz):

            for mtz_file in os.listdir(os.path.join(data_folder, xtal)):
                if "refine" in mtz_file and ".mtz" in mtz_file:
                    input_mtz = os.path.join(data_folder, xtal, mtz_file)
                else:
                    print("mtz doesn't exist, skipping: {}".format(input_pdb))
                    continue

        out_root = (
            "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/150919"
        )
        refinement_script_dir = (
            "/dls/science/groups/i04-1/elliot-dev/Work/"
            "exhaustive_search_data/191108/scripts"
        )

        # Refmac
        out_dir = os.path.join(out_root, "refmac", xtal)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_mtz = os.path.join(out_dir, "refine.mtz")
        out_pdb = os.path.join(out_dir, "refine.pdb")

        print("refmac")

        if not os.path.exists(out_pdb) and not os.path.exists(out_mtz):

            shutil.copy(src=input_pdb, dst=os.path.join(out_dir, "input.pdb"))
            shutil.copy(src=input_mtz, dst=os.path.join(out_dir, "input.mtz"))
            os.system("cd {}; giant.split_conformations input.pdb".format(out_dir))
            input_split_pdb = os.path.join(out_dir, "input.split.bound-state.pdb")

            write_refmac_csh(
                pdb=input_split_pdb,
                crystal=xtal,
                cif=input_cif,
                mtz=input_mtz,
                out_dir=out_dir,
                refinement_script_dir=refinement_script_dir,
                refinement_type="refmac",
                script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                ncyc=50,
                ccp4_path="/dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh",
            )
            csh_file = os.path.join(
                refinement_script_dir, "{}_{}.csh".format(xtal, "refmac")
            )
            os.system("qsub {}".format(csh_file))

        # Phenix
        out_dir = os.path.join(out_root, "phenix", xtal)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_mtz = os.path.join(out_dir, "refine_001.mtz")
        out_pdb = os.path.join(out_dir, "refine_001.pdb")

        print("phenix")

        print(out_pdb)
        print(os.path.exists(out_pdb))

        if not os.path.exists(out_pdb) and not os.path.exists(out_mtz):

            shutil.copy(src=input_pdb, dst=os.path.join(out_dir, "input.pdb"))
            shutil.copy(src=input_mtz, dst=os.path.join(out_dir, "input.mtz"))
            os.system("cd {}; giant.split_conformations input.pdb".format(out_dir))
            input_split_pdb = os.path.join(out_dir, "input.split.bound-state.pdb")

            if compound in issue_compounds:
                cmds = "module load phenix\n"
                cmds += "source /dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n"
                cmds += "mkdir {out_dir};\n".format(out_dir=out_dir)
                cmds += "cd {}\n".format(out_dir)
                cmds += "mkdir ./ready_set;\n"
                cmds += "cd ready_set;\n"
                cmds +=  "phenix.ready_set {pdb}\n".format(out_dir=out_dir, pdb=input_pdb)
                cmds += "cd {}\n".format(out_dir)
                input_cif = os.path.join(out_dir,
                                         "ready_set",
                                         "{pdb}.ligands.cif".format(pdb=os.path.basename(input_pdb).strip(".pdb")))
                os.system(cmds)

            write_phenix_csh(
                pdb=input_split_pdb,
                mtz=input_mtz,
                cif=input_cif,
                script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                refinement_script_dir=refinement_script_dir,
                out_dir=out_dir,
                refinement_type="bound",
                crystal=xtal,
                ncyc=20,
            )
            csh_file = os.path.join(
                refinement_script_dir, "{}_{}_bound.csh".format(xtal, "phenix")
            )
            print(csh_file)
            os.system("qsub {}".format(csh_file))

        continue

        # Buster
        out_dir = os.path.join(out_root, "buster", xtal)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_mtz = os.path.join(out_dir, "refine.mtz")
        out_pdb = os.path.join(out_dir, "refine.pdb")

        print("buster")

        if not os.path.exists(out_pdb) and not os.path.exists(out_mtz):

            shutil.copy(src=input_pdb, dst=os.path.join(out_dir, "input.pdb"))
            shutil.copy(src=input_mtz, dst=os.path.join(out_dir, "input.mtz"))
            os.system("cd {}; giant.split_conformations input.pdb".format(out_dir))
            input_split_pdb = os.path.join(out_dir, "input.split.bound-state.pdb")

            write_buster_csh(
                pdb=input_split_pdb,
                mtz=input_mtz,
                cif=input_cif,
                out_dir=out_dir,
                script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                refinement_script_dir=refinement_script_dir,
                refinement_type="bound",
                crystal=xtal,
            )
            csh_file = os.path.join(
                refinement_script_dir, "{}_{}_{}.csh".format(xtal, "buster", "bound")
            )
            os.system("qsub {}".format(csh_file))

        # exhaustive
        out_dir = os.path.join(out_root, "exhaustive", xtal)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_csv = os.path.join(out_dir, "exhaustive_search.csv")

        if not os.path.exists(out_csv):
            write_exhaustive_csh(
                pdb=input_pdb,
                mtz=input_mtz,
                script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                refinement_script_dir=refinement_script_dir,
                out_dir=out_dir,
                crystal=xtal,
                exhaustive_multiple_sampling="/dls/science/groups/i04-1/elliot-dev/Work/"
                "exhaustive_search/exhaustive/"
                "exhaustive_multiple_sampling.py",
                ccp4_path="/dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh",
            )
            csh_file = os.path.join(
                refinement_script_dir, "{}_{}.csh".format(xtal, "exhaustive")
            )

            os.system("qsub {}".format(csh_file))

        # refmac superposed
        out_dir = os.path.join(out_root, "refmac_superposed", xtal)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_mtz = os.path.join(out_dir, "refine.mtz")
        out_pdb = os.path.join(out_dir, "refine.pdb")

        if not os.path.exists(out_pdb) and not os.path.exists(out_mtz):
            shutil.copy(src=input_pdb, dst=os.path.join(out_dir, "input.pdb"))
            shutil.copy(src=input_mtz, dst=os.path.join(out_dir, "input.mtz"))

            os.system(
                "source /dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh;"
                + "cd {};".format(out_dir)
                + "giant.make_restraints {}\n".format(
                    os.path.join(out_dir, "input.pdb")
                )
            )

            with open(
                os.path.join(out_dir, "multi-state-restraints.refmac.params"), "a+"
            ) as param_file:
                param_file.write("NCYC 50")

            cmds = "source /dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n"
            cmds += "cd {}\n".format(out_dir)

            cmds += "giant.quick_refine {} {} {} params={} program={}\n".format(
                input_pdb,
                input_mtz,
                input_cif,
                "multi-state-restraints.refmac.params",
                "refmac",
            )
            cmds += "giant.split_conformations refine.pdb"

            # File location and name
            csh_file = os.path.join(
                refinement_script_dir, "{}_refmac_superposed.csh".format(xtal)
            )

            # Write file
            with open(csh_file, "w") as csh_f:
                csh_f.write(cmds)

            os.system("qsub {}".format(csh_file))


        # phenix superposed
        print("---------------------------ps---------------------------------------")
        out_dir = os.path.join(out_root, "phenix_superposed", xtal)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_mtz = os.path.join(out_dir, "refine.mtz")
        out_pdb = os.path.join(out_dir, "refine.pdb")

        if not os.path.exists(out_pdb) and not os.path.exists(out_mtz):
            shutil.copy(src=input_pdb, dst=os.path.join(out_dir, "input.pdb"))
            shutil.copy(src=input_mtz, dst=os.path.join(out_dir, "input.mtz"))

            os.system(
                "source /dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh;"
                + "cd {};".format(out_dir)
                + "giant.make_restraints {}\n".format(
                    os.path.join(out_dir, "input.pdb")
                )
            )

            with open(
                os.path.join(out_dir, "multi-state-restraints.phenix.params"), "a+"
            ) as param_file:
                param_file.write("refinement.main.number_of_macro_cycles=20")

            cmds = "module load phenix\n"
            cmds += "source /dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n"
            cmds += "mkdir {out_dir};\n".format(out_dir=out_dir)
            cmds += "cd {}\n".format(out_dir)

            if compound in issue_compounds:
                cmds += "mkdir ./ready_set;\n"
                cmds += "cd ready_set;\n"
                cmds +=  "phenix.ready_set {pdb}\n".format(out_dir=out_dir, pdb=input_pdb)
                cmds += "cd {}\n".format(out_dir)
                input_cif = os.path.join(out_dir,
                                         "ready_set",
                                         "{pdb}.ligands.cif".format(pdb=os.path.basename(input_pdb).strip(".pdb")))


            cmds += "giant.quick_refine {} {} {} params={} program={}\n".format(
                input_pdb,
                input_mtz,
                input_cif,
                os.path.join(out_dir, "multi-state-restraints.phenix.params"),
                "phenix",
            )
            cmds += "giant.split_conformations refine.pdb"

            # File location and name
            csh_file = os.path.join(
                refinement_script_dir, "{}_phenix_superposed.csh".format(xtal)
            )

            # Write file
            with open(csh_file, "w") as csh_f:
                csh_f.write(cmds)

            os.system("qsub {}".format(csh_file))


        # buster superposed
        out_dir = os.path.join(out_root, "buster_superposed", xtal)

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        out_mtz = os.path.join(out_dir, "refine.mtz")
        out_pdb = os.path.join(out_dir, "refine.pdb")

        if not os.path.exists(out_pdb) and not os.path.exists(out_mtz):
            shutil.copy(src=input_pdb, dst=os.path.join(out_dir, "input.pdb"))
            shutil.copy(src=input_mtz, dst=os.path.join(out_dir, "input.mtz"))

            cmds = "source /dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n"
            cmds += "cd {}\n".format(out_dir)
            cmds += "giant.make_restraints {}\n".format(
                os.path.join(out_dir, "input.pdb")
            )
            cmds += "giant.quick_refine {} {} {} params={} program={}\n".format(
                input_pdb, input_mtz, input_cif, "params.gelly", "buster"
            )
            cmds += "giant.split_conformations refine.pdb"

            # File location and name
            csh_file = os.path.join(
                refinement_script_dir, "{}_buster_superposed.csh".format(xtal)
            )

            # Write file
            with open(csh_file, "w") as csh_f:
                csh_f.write(cmds)

            os.system("qsub {}".format(csh_file))

