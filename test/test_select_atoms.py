import os
import unittest
from utils.phil import master_phil
from utils.select_atoms import get_bound_ground_states

import scitbx.array_family.flex

class TestSelectAtoms(unittest.TestCase):

    def setUp(self):
        """Provide setup for the test on """

        self.params = master_phil.extract()

        self.params.input.xtal_name = "FALZA-x0085"

        self.params.input.in_path = os.path.join(os.path.realpath(
            "./test/resources"), self.params.input.xtal_name)

        self.params.input.pdb = os.path.join(self.params.input.in_path, "refine.pdb")

    def test_get_bound_ground_states(self):

        ground_states, bound_states = get_bound_ground_states(pdb=self.params.input.pdb,
                                                              params=self.params)
        #Type checking
        self.assertIsInstance(ground_states, list)
        self.assertIsInstance(ground_states[0], list)
        self.assertIsInstance(ground_states[0][1], int)
        self.assertIsInstance(ground_states[0][0],
                              scitbx.array_family.flex.bool)

        self.assertIsInstance(bound_states, list)
        self.assertIsInstance(bound_states[0], list)
        self.assertIsInstance(bound_states[0][1], int)
        self.assertIsInstance(bound_states[0][0],
                              scitbx.array_family.flex.bool)

        # Check for specific pdb file parsing
        self.assertEqual(ground_states[0][1], 1)
        self.assertEqual(len(ground_states), 12)

        self.assertEqual(bound_states[0][1], 1)
        self.assertEqual(len(bound_states), 10)





if __name__ == '__main__':
    unittest.main()
