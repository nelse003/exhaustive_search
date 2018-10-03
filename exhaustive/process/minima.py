from __future__ import division, print_function

import logging
import pandas as pd
import iotbx
from iotbx.pdb import hierarchy

from exhaustive.utils.utils import process_refined_pdb_bound_ground_states
from exhaustive.utils.utils import u_iso_to_b_fac, b_fac_to_u_iso, get_minimum_fofc

logging = logging.getLogger(__name__)

def get_minima_within_b_range(csv_file, b_lower, b_upper):

    df = pd.read_csv(csv_file)
    u_lower = b_fac_to_u_iso(b_lower)
    u_upper = b_fac_to_u_iso(b_upper)
    print(df)


def write_minima_pdb(input_pdb,output_pdb,csv_name, params):

    min_occ, min_u_iso, _ = get_minimum_fofc(csv_name)

    bound_states, ground_states = process_refined_pdb_bound_ground_states(
        input_pdb, params)
    pdb_inp = iotbx.pdb.input(input_pdb)
    hier = pdb_inp.construct_hierarchy()

    for chain in hier.only_model().chains():
        for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                for atom in atom_group.atoms():

                    for ground_state in ground_states:
                        num_altlocs = ground_state[1]
                        if ground_state[0][atom.i_seq]:
                            atom.occ = (1 - min_occ) / num_altlocs
                            atom.b = u_iso_to_b_fac(min_u_iso)

                    for bound_state in bound_states:
                        num_altlocs = bound_state[1]
                        if bound_state[0][atom.i_seq]:
                            atom.set_occ(min_occ/num_altlocs)
                            atom.set_b(u_iso_to_b_fac(min_u_iso))


    with open(output_pdb,"w") as f:
        f.write(hier.as_pdb_string(
            crystal_symmetry=hierarchy.input(input_pdb).crystal_symmetry()))


get_minima_within_b_range(csv_file="/dls/science/groups/i04-1/elliot-dev/Work/"\
                          "exhaustive_search_data/covalent_ratios_exhaus_sep_18/NUDT7A-x6208",
                          b_lower=0.3,
                          b_upper=0.7)
