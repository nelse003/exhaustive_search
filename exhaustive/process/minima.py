import csv
import os
import sys
import logging

import iotbx
import libtbx.phil
import numpy as np
from iotbx.pdb import hierarchy

from exhaustive.utils.select import process_refined_pdb_bound_ground_states
from exhaustive.utils.utils import u_iso_to_b_fac
from exhaustive.validation.Repeating_exhaustive_search import get_in_refinement_or_better

#################################################
master_phil = libtbx.phil.parse("""
input{
    database_path = "/dls/labxchem/data/2018/lb18145-55/processing/database/soakDBDataFile.sqlite"
        .type = path
    csv_name = 'u_iso_occupancy_vary_new_atoms'
        .type = str
}
output{
    out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/occupancy_group_with_refinement"
        .type = str
    minima_csv_name = "min_occ_u_iso_NUDT22_with_refinement"
        .type = str
}
options{
    cat = "cat"
        .type = str
}
""", process_includes=True)
###################################################

logger = logging.getLogger(__name__)

def write_minima_pdb(input_pdb,output_pdb,csv_name, params):

    min_occ, min_u_iso, _ = get_minimum_fofc(csv_name)

    bound_states, ground_states = process_refined_pdb_bound_ground_states(input_pdb, params)
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


    with open(output_pdb,"w") as file:
        file.write(hier.as_pdb_string(crystal_symmetry=hierarchy.input(input_pdb).crystal_symmetry()))



def run(params):

    pass

    # print(params.input.database_path)
    # get_all_minima(params)
    # print("AAAAAAAAAAAAA")
    # input_pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1787/refine.pdb"
    # csv_name = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/NUDT7_Copied_atoms/NUDT7A-x1787/u_iso_occupancy_vary_new_atoms"
    # write_minima_pdb(input_pdb, csv_name)


if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:])