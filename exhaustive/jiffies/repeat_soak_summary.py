import os
import sys
import pandas as pd

sys.path.append("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search")
from exhaustive.utils.utils_ccp4 import read_ligand
from exhaustive.utils.utils import get_minimum_fofc, u_iso_to_b_fac

if __name__ == "__main__":

    """Summarise occupancy stats from repeat soak refienements
    
    Uses ccp4-python
    """

    out_root = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/150919"

    xtals_with_compound = {}

    # DCP2B FMOPL00213a
    # What are these hits?

    #DCP2B FMOPL000435a
    xtals_with_compound["DCP2B-x0146"] = 'FMOPL000435a'
    start_xtal_num = 993
    end_xtal_num = 1048
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "DCP2B-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'FMOPL000435a'

    # NUDT7A OX210 1774-1810, 0299
    xtals_with_compound['NUDT7A-x0299'] = 'OX210'
    start_xtal_num = 1774
    end_xtal_num = 1810
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT7A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'OX210'

    # These are missing in refmac...

    #NUDT7A NUOOA000181a
    #1739 - 1773
    #2006 - 2073
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


    # These are missing in refmac...

    # NUDT22A FMOPL00622a
    #from DSI poised library (x0938-x0976)
    xtals_with_compound['NUDT22A-x0182'] = 'FMOPL000622a_DSPL'
    start_xtal_num = 909
    end_xtal_num = 937
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] = 'FMOPL000622a_DSPL'

    # NUDT22A FMOPL00622a
    #from DSI poised library (x0938-x0976)
    start_xtal_num = 938
    end_xtal_num = 976
    for num in range(start_xtal_num, end_xtal_num + 1):
        xtal_name = "NUDT22A-x" + "{0:0>4}".format(num)
        xtals_with_compound[xtal_name] =' FMOPL000622a_DSI_poised'

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

    occ_df_list =[]
    for xtal, compound in xtals_with_compound.items():

        for program_folder in os.listdir(out_root):

            if not os.path.isdir(os.path.join(out_root, program_folder)):
                continue
            elif program_folder=="scripts":
                continue

            if "phenix" in program_folder and not "superposed" in program_folder:
                pdb = os.path.join(out_root, program_folder, xtal,"refine_001.pdb")

            elif program_folder=="buster" or program_folder=="refmac":
                pdb = os.path.join(out_root, program_folder, xtal, "refine.pdb")

            # update split confs
            elif program_folder=="phenix_superposed":
                pdb = os.path.join(out_root, program_folder, xtal, "refine.pdb")
            else:
                pdb = os.path.join(out_root, program_folder, xtal, "refine.split.bound-state.pdb")

            if not "exhaustive" in program_folder:
                if os.path.exists(pdb):
                    occ_df = read_ligand(pdb_path=pdb)
                else:
                    continue
            else:
                csv = os.path.join(out_root,
                                    program_folder,
                                    xtal,
                                    "exhaustive_search.csv")

                if not os.path.exists(csv):
                    continue

                occ, u_iso, fo_fc = get_minimum_fofc(csv_name=csv)
                exhaustive_data = {"Occupancy": occ,
                        "B_factor":u_iso_to_b_fac(u_iso),
                        "min_fo_fc": fo_fc}

                occ_df = pd.DataFrame(data=exhaustive_data, index=[0])


                # Dealin with two copies of ligand
                if compound == "OX210":
                    occ_df_2 = occ_df.copy(deep=True)
                    occ_df['resid'] = 1
                    occ_df_2['resid'] = 2
                    occ_df = pd.concat([occ_df, occ_df_2])
                else:
                    occ_df['resid'] = 1


            occ_df['xtal'] = xtal
            occ_df['refinement'] = program_folder
            occ_df['compound'] = compound

            occ_df_list.append(occ_df)

            if xtal=="NUDT22A-x1038" and program_folder == "phenix_superposed":
                print(pdb)
                print(occ_df)


    occ_all_df = pd.concat(occ_df_list)
    occ_all_df.to_csv("/dls/science/groups/i04-1/elliot-dev/Work/"
                        "repeat_soak_plots/150919/occ_all.csv")


