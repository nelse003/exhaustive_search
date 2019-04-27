import os
import unittest
import shutil

from exhaustive import run as exhaustive
from utils.phil import master_phil
from utils.utils import get_minimum_fofc


# TODO Write TestComputeMaps
class TestComputeMaps(unittest.TestCase):
    """
    Test the ability to compute maps using compute_maps()
    """

    def setUp(self):
        self.params = master_phil.extract()

    @unittest.skip("Not implemented")
    def test_compute_maps(self):
        compute_maps()


class TestExhaustiveSearch(unittest.TestCase):
    """
    Test the main loop of exhaustive search.
    """

    def setUp(self):
        """Provide the minimal number of parameters to test exhaustive search"""

        self.params = master_phil.extract()

        self.params.input.xtal_name = "FALZA-x0085"
        self.params.input.in_path = os.path.join(
            os.path.realpath("./test/resources"), self.params.input.xtal_name
        )
        self.params.validate.input.base_mtz = os.path.join(
            self.params.input.in_path, "FALZA-x0085.free.mtz"
        )
        self.params.input.mtz = os.path.join(
            self.params.input.in_path, "FALZA-x0085.free.mtz"
        )
        self.params.input.pdb = os.path.join(self.params.input.in_path, "refine.pdb")
        self.params.output.out_dir = os.path.realpath("./test/output")
        self.params.output.log_dir = os.path.realpath(
            os.path.join("./test/output", "logs")
        )
        self.params.exhaustive.options.step = 0.05

    def tearDown(self):
        """Remove test files"""
        # TODO Fix Teardown method

        # This Fails due to a nfs error with files not being closed.

        # OSError: [Errno 16] Device or resource busy:
        # '/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/test/output/logs/.nfs00000000a0d05de80000319f'

        shutil.rmtree(os.path.realpath("./test/output"))

    def test_exhaustive_search(self):
        """ Test with minimal number of parameters changed from default."""

        self.params.exhaustive.output.csv_name = os.path.join(
            self.params.output.out_dir, "test.csv"
        )
        self.params.exhaustive.options.step = 0.2
        exhaustive(self.params)
        bound_occ, u_iso, fofc = get_minimum_fofc(
            self.params.exhaustive.output.csv_name
        )
        self.assertAlmostEqual(0.6, bound_occ)
        self.assertAlmostEqual(0.4, u_iso)

    @unittest.skip("Freezing")
    def test_exhaustive_search_multiprocessing(self):
        """ Test with minimal number of parameters changed from default."""

        self.params.settings.processes = 4

        self.params.exhaustive.output.csv_name = os.path.join(
            self.params.output.out_dir, "test.csv"
        )
        self.params.exhaustive.options.step = 0.2
        exhaustive(self.params)
        bound_occ, u_iso, fofc = get_minimum_fofc(
            self.params.exhaustive.output.csv_name
        )
        self.assertAlmostEqual(0.6, bound_occ)
        self.assertAlmostEqual(0.4, u_iso)

    def test_convex_hull_exhaustive_search(self):
        self.params.exhaustive.output.csv_name = os.path.join(
            self.params.output.out_dir, "test.csv"
        )
        exhaustive(self.params)
        bound_occ, u_iso, fofc = get_minimum_fofc(
            self.params.exhaustive.output.csv_name
        )
        self.assertAlmostEqual(0.6, bound_occ)
        self.assertAlmostEqual(0.35, u_iso)

    @unittest.skip("Not implemented")
    def test_convex_hull_buffer_exhaustive_search(self):
        self.assertEqual(True, False)

    @unittest.skip("Not implemented")
    def test_convex_hull_buffer_ignore_nearest(self):
        self.assertEqual(True, False)

    @unittest.skip("Not implemented")
    def test_step_size(self):
        self.assertEqual(True, False)


# TODO Write TestOccupancyCartPoints
class TestOccupancyCartPoints(unittest.TestCase):
    """
    Test the generation of cartesian points to loop over

    What mock objects does this need?
    """

    @unittest.skip("Not implemented")
    def test_occupancy_cart_points(self):
        self.assertEqual(True, False)


# TODO Write TestCalculateMeanFofcCartSites
class TestCalculateMeanFofcCartSites(unittest.TestCase):
    """
    Test the calculation of |Fo-Fc| given cartesian sitess
    """

    @unittest.skip("Not implemented")
    def test_calc_mean_fofc_cart_sites(self):
        self.assertEqual(True, False)


if __name__ == "__main__":
    unittest.main()
