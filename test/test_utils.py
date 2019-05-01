import unittest
from exhaustive.utils.utils import b_to_u_iso
from exhaustive.utils.utils import u_iso_to_b_fac


class TestB_fac(unittest.TestCase):
    def test_b_to_u_iso(self):
        """Test the calculation of u_iso from B factor"""
        u_iso = b_to_u_iso(b_fac=40)
        self.assertAlmostEqual(u_iso, 0.5066059182116889, places=5)


class TestU_iso(unittest.TestCase):
    def test_b_to_u_iso(self):
        """Test the calculation of B factor from u_iso"""
        b_fac = u_iso_to_b_fac(u_iso=0.5)
        self.assertAlmostEqual(b_fac, 39.47841760435743, places=5)


if __name__ == "__main__":
    unittest.main()
