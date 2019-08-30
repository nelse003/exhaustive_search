import sys
# needed to call from ccp4 python
sys.path.append("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search")

import argparse
import iotbx
from iotbx.pdb import hierarchy
from exhaustive.utils.select_atoms import get_bound_ground_states
from exhaustive.utils.utils import get_minimum_fofc
from exhaustive.utils.utils import u_iso_to_b_fac
from phil import master_phil


def write_minima_pdb(input_pdb, output_pdb, csv_name, params):
    """
    Write pdb from the minima in exhaustive search

    Parameters
    ----------
    input_pdb: str
        path to input pdb to take structure from

    output_pdb: str
        path to write strucutre to

    csv_name: str
        path to exhaustive search csv

    params: str
        parameter

    Returns
    -------

    """

    min_occ, min_u_iso, _ = get_minimum_fofc(csv_name)

    bound_states, ground_states = get_bound_ground_states(input_pdb, params)
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
                            atom.set_occ(min_occ / num_altlocs)
                            atom.set_b(u_iso_to_b_fac(min_u_iso))

    with open(output_pdb, "w") as f:
        f.write(
            hier.as_pdb_string(
                crystal_symmetry=hierarchy.input(input_pdb).crystal_symmetry()
            )
        )

if __name__ == "__main__":

    params = master_phil.extract()

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_pdb")
    parser.add_argument("-o", "--output_pdb")
    parser.add_argument("-c", "--csv_name")

    args = parser.parse_args()

    write_minima_pdb(input_pdb=args.input_pdb,
                     output_pdb=args.output_pdb,
                     csv_name=args.csv_name,
                     params=params)