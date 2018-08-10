from utils.utils import remove_residues

remove_residues("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"
             "covalent_ratios_phenix/NUDT7A-x1906/refine.pdb",
             "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"
             "covalent_ratios_phenix/NUDT7A-x1906/refine_removed.pdb",
                [['73','A','C'],['73','A','D'],['1','E','D']])