import os
import unittest
import cctbx
import mmtbx

from xtal_model_data import XtalModelData
from utils.phil import master_phil

class TestXtalModelData(unittest.TestCase):

    def setUp(self):
        """Provide setup for the test on """

        self.params = master_phil.extract()

        self.params.input.xtal_name = "FALZA-x0085"

        self.params.input.in_path = os.path.join(os.path.realpath(
            "./test/resources"), self.params.input.xtal_name)

        self.params.validate.input.base_mtz = os.path.join(self.params.input.in_path,
                                                           "FALZA-x0085.free.mtz")

        self.params.input.mtz = os.path.join(self.params.input.in_path,
                                             "FALZA-x0085.free.mtz")

        self.params.input.pdb = os.path.join(self.params.input.in_path, "refine.pdb")

    def test_init(self):
        """Check the class constructor

        Checks that objects of the correct type are produced.
        """
        xtal_model_data = XtalModelData(self.params)

        assert xtal_model_data.pdb == self.params.input.pdb

        assert xtal_model_data.mtz == self.params.input.mtz

        self.assertIsInstance(xtal_model_data.xrs,
                              cctbx.xray.structure)

        self.assertIsInstance(xtal_model_data.inputs,
                              mmtbx.utils.process_command_line_args)

        self.assertIsInstance(xtal_model_data.crystal_gridding,
                              cctbx.maptbx.crystal_gridding)

        self.assertIsInstance(xtal_model_data.fmodel,
                              mmtbx.f_model.f_model.manager)

    def test_determine_states(self):

        xtal_model_data = XtalModelData(self.params)
        self.assertIsInstance(xtal_model_data.ground_states,
                              list)

        assert True == False

    def test_iter_u_iso_occ(self):

        xtal_model_data = XtalModelData(self.params)
        u_iso_occ = xtal_model_data.iter_u_iso_occ()

        self.assertIsInstance(u_iso_occ, list)

        assert True == False

    def test_calculate_mean_fofc(self):

        assert True == False