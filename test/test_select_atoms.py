import os
import unittest

import iotbx
import scitbx.array_family.flex

from utils.phil import master_phil
from utils.select_atoms import get_bound_ground_states
from utils.select_atoms import get_selection_altloc_resseq_chain


class TestSelectAtoms(unittest.TestCase):

    def setUp(self):
        """Provide setup for the test on """

        self.params = master_phil.extract()

        self.params.input.xtal_name = "FALZA-x0085"

        self.params.input.in_path = os.path.join(os.path.realpath(
            "./test/resources"), self.params.input.xtal_name)

        self.params.input.pdb = os.path.join(self.params.input.in_path, "refine.pdb")

        pdb_inp = iotbx.pdb.input(self.params.input.pdb )
        hier = pdb_inp.construct_hierarchy()
        self.sel_cache = hier.atom_selection_cache()

    def test_get_bound_ground_states(self):

        ground_states, bound_states = get_bound_ground_states(pdb=self.params.input.pdb,
                                                              params=self.params)
        # Type checking
        self.assertIsInstance(ground_states, list)
        self.assertIsInstance(ground_states[0], tuple)
        self.assertIsInstance(ground_states[0][1], int)
        self.assertIsInstance(ground_states[0][0],
                              scitbx.array_family.flex.bool)

        self.assertIsInstance(bound_states, list)
        self.assertIsInstance(bound_states[0], tuple)
        self.assertIsInstance(bound_states[0][1], int)
        self.assertIsInstance(bound_states[0][0],
                              scitbx.array_family.flex.bool)

        # Check for specific pdb file parsing
        self.assertEqual(ground_states[0][1], 2)
        self.assertEqual(len(ground_states), 6)

        self.assertEqual(bound_states[0][1], 2)
        self.assertEqual(len(bound_states), 5)


    def test_get_selection_altloc_resseq_chain(self):

        sel_alt = get_selection_altloc_resseq_chain(sel_cache=self.sel_cache,
                                          altlocs=['A','B'],
                                          resseq=120,
                                          chain='S')

        # Type checking
        self.assertIsInstance(sel_alt, tuple)
        self.assertIsInstance(sel_alt[0], scitbx.array_family.flex.bool)
        self.assertIsInstance(sel_alt[1],int)

        # Check for specific pdb file parsing
        self.assertEqual(sel_alt[1], 2)

if __name__ == '__main__':
    unittest.main()
