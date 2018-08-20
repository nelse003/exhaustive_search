import unittest
from exhaustive.phil import master_phil
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
    
    What mock objects does this need?
    """
    def setUp(self):

        """Provide a """

        self.params = master_phil.extract()
        self.params.input.pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/" \
                                "initial_model/NUDT7A-x0299/refine.pdb"
        self.params.input.mtz = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/" \
                                "initial_model/NUDT7A-x0299/refine.mtz"
        self.params.input.xtal_name = "NUDT7A-x0299"
        self.params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                                     "exhaustive_search_data/tests/"
        self.params.output.log_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                                     "exhaustive_search_data/tests/logs"

    def test_exhaustive_search(self):

        self.params.exhaustive.output.csv_name = os.path.join(self.params.output.out_dir, "test.csv")
        exhaustive(self.params)
        bound_occ, u_iso, fofc = get_minimum_fofc(self.params.exhaustive.output.csv_name)
        self.assertAlmostEqual(0.7,bound_occ)
        self.assertAlmostEqual(0.45,u_iso)

    def test_convex_hull_exhaustive_search(self):
        self.assertEqual(True, False)

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
