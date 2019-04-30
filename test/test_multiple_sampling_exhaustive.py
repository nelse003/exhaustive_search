import os

from exhaustive.exhaustive_multiple_sampling import run as multiple_exhaustive
from test_exhaustive_search import TestExhaustiveSearch
from exhaustive.utils.utils import get_minimum_fofc


class TestMultipleSamplingExhaustiveSearch(TestExhaustiveSearch):
    """
    Test the multiple sampling exhaustive search method.

    Inherits setup method from TestExhaustiveSearch
    """

    def test_multiple_exhaustive_search(self):
        """ Test with minimal number of parameters changed from default."""

        self.params.exhaustive.output.csv_name = os.path.join(
            self.params.output.out_dir, "test.csv"
        )
        multiple_exhaustive(self.params)
        bound_occ, u_iso, fofc = get_minimum_fofc(
            self.params.exhaustive.output.csv_name
        )
        self.assertAlmostEqual(0.6, bound_occ)
        self.assertAlmostEqual(0.33, u_iso)
