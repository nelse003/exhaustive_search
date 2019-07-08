import os
import sys

sys.path.append("/dls/science/groups/i04-1/elliot-dev/parse_xchemdb")
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

    # DCP2B FMOPL000435a
    xtals_with_compound["DCP2BA-x0146"] = 'FMOPL000435a'
    start_xtal_num = 993
    end_xtal_num = 1048
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "DCP2BA-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'FMOPL000435a'

    # NUDT7A OX210 1774-1810, 0299
    xtals_with_compound['NUDT7A-x0299'] = 'OX210'
    start_xtal_num = 1774
    end_xtal_num = 1810
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT7A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'OX210'

    # NUDT7A NUOOA000181a
    # 1739 - 1773
    # 2006 - 2073
    start_xtal_num = 1739
    end_xtal_num = 1773
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT7A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'NUOOA000181a'
    start_xtal_num = 2006
    end_xtal_num = 2073
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT7A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'NUOOA000181a'


    # # NUDT22A FMOPL00622a
    # from DSI poised library (x0938-x0976)
    xtals_with_compound['NUDT22A-x0182'] = 'FMOPL000622a_DSPL'
    start_xtal_num = 909
    end_xtal_num = 937
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'FMOPL000622a_DSPL'

    # # NUDT22A FMOPL00622a
    # from DSI poised library (x0938-x0976)
    start_xtal_num = 938
    end_xtal_num = 976
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'FMOPL000622a_DSI_poised'

    # NUDT22A 133725a x0421, x1040 - x1059
    start_xtal_num = 1040
    end_xtal_num = 1059
    xtals_with_compound['NUDT22A-x0421'] = 'N13725a'
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'N13725a'

    # NUDT22A 13369a x0243, x977 to x1008

    start_xtal_num = 977
    end_xtal_num = 1008
    xtals_with_compound['NUDT22A-x0243'] = 'N13369a'
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'N13369a'

    # NUDT22A 13663a x0391, x1009 to x1039
    xtals_with_compound['NUDT22A-x0391'] = 'N13663a'
    start_xtal_num = 1009
    end_xtal_num = 1039
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'N13663a'

    compound_codes = set()
    for compounds in xtals_with_compound.values():
        compound_codes.add(compounds)

    # Define the file name and the folders to loop over
    folders = {
        "phenix/": ("refine_001.pdb", "refine_001.mtz", "phenix"),
        "refmac/": ("refine.pdb", "refine.mtz", "refmac"),
        "buster/": ("refine.pdb", "refine.mtz", "buster"),
        "buster_superposed": ("refine.pdb", "refine.mtz", "buster_superposed"),
        "exhaustive" :()
    }

    for xtal, compound in xtals_with_compound.items():

        if "DCP2B" in xtal:
            data_folder =  "/dls/science/groups/i04-1/elliot-dev/Work/" \
                           "exhaustive_search_data/DCP2B_18_09_20_exhaus"
            input_cif =

        elif "NUDT22" in xtal:
            data_folder_pre = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                          "exhaustive_search_data/NUDT22_repeat_soaks"

            if compound == "FMOPL000622a_DSI_poised":
                folder_ext = "FMOPL00622a: Alternate Diastereomer"
            elif compound == "FMOPL000622a_DSPL":
                folder_ext = "FMOPL000622a"
            else:
                folder_ext = compound

            data_folder = os.path.join(data_folder_pre, folder_ext)

        elif "NUDT7" in xtal:
            if compound == "OX210":
                data_folder = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                              "exhaustive_search_data/NUDT7_Copied_atoms/OX210"

                input_cif = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/" \
                            "NUDT7_Copied_atoms/OX210/NUDT7A-x1787/OX-210.cif"
            else:
                # This isn't copied atoms. Can't find that...
                data_folder ="/dls/labxchem/data/2017/lb18145-3/processing/analysis/initial_model"


        input_pdb = os.path.join(data_folder, xtal, "refine.pdb")
        input_mtz = os.path.join(data_folder, xtal, "refine.mtz")

        if os.path.exists(input_pdb):
            pass
            #print("PDB exists: {}".format(input_pdb))
        else:
            print("PDB doesn't exist: {}".format(input_pdb))

        continue
        exit()
        for folder in folders:
            for crystal_folder in os.listdir(os.path.join(root, folder)):
                if os.path.isdir(os.path.join(root, folder, crystal_folder)):

                    exit()

                    if args.program == "refmac":
                        write_refmac_csh(
                            pdb=pdb_adj,
                            crystal=xtal,
                            cif=cif,
                            mtz=mtz,
                            out_dir=os.path.join(out_dir, xtal),
                            refinement_script_dir=refinement_script_dir,
                            script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                            ncyc=50,
                            ccp4_path="/dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh",
                        )
                        csh_file = os.path.join(
                            refinement_script_dir, "{}_{}.csh".format(xtal, "bound")
                        )

                    elif args.program == "phenix":
                        write_phenix_csh(
                            pdb=pdb_adj,
                            mtz=mtz,
                            cif=cif,
                            script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                            refinement_script_dir=refinement_script_dir,
                            out_dir=os.path.join(out_dir, xtal),
                            crystal=xtal,
                            ncyc=20,
                        )
                        csh_file = os.path.join(
                            refinement_script_dir, "{}_{}.csh".format(xtal, "phenix")
                        )

                    elif args.program == "buster":
                        write_buster_csh(
                            pdb=pdb_adj,
                            mtz=mtz,
                            cif=cif,
                            out_dir=os.path.join(out_dir, xtal),
                            script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                            refinement_script_dir=refinement_script_dir,
                            crystal=xtal,
                        )
                        csh_file = os.path.join(
                            refinement_script_dir, "{}_{}.csh".format(xtal, "buster")
                        )

                    elif args.program == "buster_b_none":
                        print(out_dir)
                        write_buster_csh(
                            pdb=pdb_adj,
                            mtz=mtz,
                            cif=cif,
                            out_dir=os.path.join(out_dir, xtal),
                            script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                            refinement_script_dir=refinement_script_dir,
                            crystal=xtal,
                            template_name="buster_b_none_template.csh",
                        )
                        csh_file = os.path.join(
                            refinement_script_dir, "{}_{}.csh".format(xtal, "buster")
                        )

                    elif args.program == "exhaustive":
                        write_exhaustive_csh(
                            pdb=refine_pdb,
                            mtz=mtz,
                            script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                            refinement_script_dir=refinement_script_dir,
                            out_dir=os.path.join(out_dir, xtal),
                            crystal=xtal,
                            exhaustive_multiple_sampling="/dls/science/groups/i04-1/elliot-dev/Work/"
                                                         "exhaustive_search/exhaustive/"
                                                         "exhaustive_multiple_sampling.py",
                            ccp4_path="/dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh",
                        )
                        csh_file = os.path.join(
                            refinement_script_dir, "{}_{}.csh".format(xtal, "exhaustive")
                        )

                        #os.system("qsub {}".format(csh_file))