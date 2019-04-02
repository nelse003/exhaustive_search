import unittest
from exhaustive.exhaustive.utils.phil import master_phil
from exhaustive.exhaustive.exhaustive import run as exhaustive
from exhaustive.exhaustive.utils.utils import get_minimum_fofc
import os

# TODO Write TestComputeMaps
class TestComputeMaps(unittest.TestCase):
    """
    Test the ability to compute maps using compute_maps()
    """

    def test_compute_maps(self):
        self.assertEqual(True, False)


# TODO Write TestExhaustiveSearch
class TestExhaustiveSearch(unittest.TestCase):
    """
    Test the main loop of exhaustive search.
    """
    def setUp(self):

        """Provide the minimal number of parameters to test exhaustive search"""

        #TODO Convert to relative import of test resources

        self.params = master_phil.extract()

        self.params.input.xtal_name = "FALZA-x0085"
        self.params.input.in_path = os.path.join(os.path.realpath(
            "./exhaustive/test/resources"), self.params.input.xtal_name)
        self.params.validate.input.base_mtz = os.path.join(self.params.input.in_path,
                                        "FALZA-x0085.free.mtz")
        self.params.input.mtz = os.path.join(self.params.input.in_path,
                                        "FALZA-x0085.free.mtz")
        self.params.input.pdb = os.path.join(self.params.input.in_path,"refine.pdb")
        self.params.output.out_dir = os.path.realpath("./exhaustive/test/output")
        self.params.output.log_dir = os.path.realpath(os.path.join("./exhaustive/test/output","logs"))

    def test_exhaustive_search(self):

        """ Test with minimal number of parameters changed from default."""

        self.params.exhaustive.output.csv_name = os.path.join(self.params.output.out_dir, "test.csv")
        exhaustive(self.params)
        bound_occ, u_iso, fofc = get_minimum_fofc(self.params.exhaustive.output.csv_name)
        self.assertAlmostEqual(0.7,bound_occ)
        self.assertAlmostEqual(0.45,u_iso)

    def test_convex_hull_exhaustive_search(self):

        self.params.exhaustive.output.csv_name = os.path.join(self.params.output.out_dir, "test.csv")
        exhaustive(self.params)
        bound_occ, u_iso, fofc = get_minimum_fofc(self.params.exhaustive.output.csv_name)
        self.assertAlmostEqual(0.7,bound_occ)
        self.assertAlmostEqual(0.45,u_iso)

    def test_convex_hull_buffer_exhaustive_search(self):
        self.assertEqual(True, False)

    def test_convex_hull_buffer_ignore_nearest(self):
        self.assertEqual(True, False)

    def test_step_size(self):
        self.assertEqual(True, False)

    def test_multiprocess_exhaustive_search(self):
        self.assertEqual(True, False)


# TODO Write TestOccupancyCartPoints
class TestOccupancyCartPoints(unittest.TestCase):
    """
    Test the generation of cartesian points to loop over

    What mock objects does this need?
    """

    def test_occupancy_cart_points(self):
        self.assertEqual(True, False)


# TODO Write TestCalculateMeanFofcCartSites
class TestCalculateMeanFofcCartSites(unittest.TestCase):
    """
    Test the calculation of |Fo-Fc| given cartesian sitess
    """

    def test_calc_mean_fofc_cart_sites(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
