import os
import unittest
import shutil

from xtal_model_data import XtalModelData
from occ_b_loop import OccBLoopCaller
from utils.phil import master_phil


class TestOccBLoop(unittest.TestCase):
    def setUp(self):
        """Provide setup for the test on """

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
        self.params.exhaustive.options.generate_mtz = False

    def tearDown(self):
        """Remove test files"""
        # TODO Fix Teardown method on diamond

        # This Fails due to a nfs error with files not being closed.

        # OSError: [Errno 16] Device or resource busy:
        # '/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/test/output/logs/.nfs00000000a0d05de80000319f'
        if os.path.exists(os.path.realpath("./test/output")):
            shutil.rmtree(os.path.realpath("./test/output"))

    def test_call(self):
        """Test a single iteration of the call method from OccBLoopCaller"""
        self.xtal_model_data = XtalModelData(self.params)
        occ_b_loop = OccBLoopCaller(xtal_model_data=self.xtal_model_data)
        result = map(occ_b_loop, [(0.6, 0.4)])
        bound_occupancy = result[0][0]
        ground_occupancy = result[0][1]
        u_iso = result[0][2]
        mean_fofc = result[0][3]

        self.assertEqual(bound_occupancy, 0.6)
        self.assertEqual(ground_occupancy, 0.4)
        self.assertEqual(u_iso, 0.4)
        self.assertAlmostEqual(mean_fofc, 0.181, places=2)

    def test_call_mtz(self):
        """Test call of OccBLoopCall generating a mtz file"""

        self.params.exhaustive.options.generate_mtz = True
        self.xtal_model_data = XtalModelData(self.params)

        occ_b_loop = OccBLoopCaller(xtal_model_data=self.xtal_model_data)
        result = map(occ_b_loop, [(0.6, 0.4)])

        assert os.path.exists(
            os.path.join(self.params.output.out_dir, "testing_0_6_0_4.mtz")
        )

    def test_call_map(self):
        """Test call of OccBLoopCall generating a map file"""

        self.params.exhaustive.options.generate_mtz = True
        self.params.exhaustive.options.generate_map = True
        self.xtal_model_data = XtalModelData(self.params)

        occ_b_loop = OccBLoopCaller(xtal_model_data=self.xtal_model_data)
        result = map(occ_b_loop, [(0.6, 0.4)])

        assert os.path.exists(
            os.path.join(self.params.output.out_dir, "testing_0_6_0_4_mFo-DFc.ccp4")
        )


if __name__ == "__main__":
    unittest.main()
