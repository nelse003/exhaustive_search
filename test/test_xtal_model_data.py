import os
import unittest

from xtal_model_data import XtalModelData
from utils.phil import master_phil

class TestXtalModelData(unittest.TestCase):

    def setUp(self):

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
        xtal_model_data = XtalModelData(self.params)

        assert xtal_model_data.pdb == self.params.input.pdb
        assert xtal_model_data.mtz == self.params.input.mtz



