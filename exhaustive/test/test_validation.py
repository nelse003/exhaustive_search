import os

import numpy as np
import unittest

from exhaustive.exhaustive.utils.utils import get_minimum_fofc, b_to_u_iso, u_iso_to_b_fac
from exhaustive.validation.validation import run as validate
from phil import master_phil


class TestExhaustiveSearch(unittest.TestCase):
    """
    Test the main loop of exhaustive search.

    What mock objects does this need?
    """

    def setUp(self):

        self.params =  master_phil.extract()
        self.params.input.xtal_name = "FALZA-x0085"
        self.params.input.in_path = os.path.join(os.path.realpath(
            "./exhaustive/test/resources"), self.params.input.xtal_name)
        self.params.input.mtz = os.path.join(self.params.input.in_path,
                                        "FALZA-x0085.free.mtz")
        self.params.input.pdb = os.path.join(self.params.input.in_path,"refine.pdb")
        self.params.output.out_dir = os.path.realpath("./exhaustive/test/output")
        self.params.output.log_dir = os.path.join(self.params.output.out_dir, "logs")
        self.params.validate.input.ground_state_pdb_path = os.path.join(
            self.params.input.in_path, "refine.output.ground-state.pdb")
        self.params.validate.input.bound_state_pdb_path = os.path.join(
            self.params.input.in_path, "refine.output.bound-state.pdb")
        self.params.validate.options.set_b = 40
        self.params.exhaustive.options.generate_mtz = False
        self.params.validate.options.use_qsub = False
        self.params.validate.options.step_simulation = 0.1
        self.params.validate.options.overwrite = True
        self.params.exhaustive.options.step = 0.1
        self.params.settings.processes = 20

        if not os.path.exists(self.params.output.out_dir):
            os.mkdir(self.params.output.out_dir)

        if not os.path.exists(self.params.output.log_dir):
            os.mkdir(self.params.output.log_dir)


    def test_validation(self):
        validate(self.params)
        check_validate_result(self.params)
        check_validate_ouput_files(self.params)


    def check_validate_result(self):

        """ Check simulated validation occupancy for all simulated occupancies"""

        for simul_occ in np.arange(self.params.validate.options.start_simul_occ,
                              self.params.validate.options.end_simul_occ,
                              self.params.validate.options.step_simulation):
            validate_occ_check(simul_occ)

    def validate_occ_check(simul_occ):

        """ Compare simulated occupancy to occupancy minima
    
        Assert statements will give AssertionError if minima is too far
        (2* step size (default 0.02)) from
        """
        occ, u_iso, fo_fc = get_minimum_fofc(self.params.output.out_dir,
                                             self.params.exhaustive.output.csv_prefix
                                             + "_occ_{}_b_{}.csv".format(
                                                 str(simul_occ).replace(".","_"),
                                                 str(self.params.validate.options.set_b)
                                             ))
        assert simul_occ - occ >= self.params.exhaustive.options.step*2,\
            "Occupancy minima {} is too farfrom simulated occupancy {}".format(
                occ,
                simul_occ)
        assert u_iso - b_to_u_iso(self.params.validate.options.set_b) \
               >= self.params.exhaustive.options.step*2, \
            "B factor minima {} is too far " \
            "from target simulated b factor {}".format(u_iso_to_b_fac(u_iso),
                                                    self.params.validate.options.set_b)


def check_validate_ouput_files(self):

    """
    Check whether file output matches expected file output.
    All files expected exist, and no extra files generated

    :param params:
    :return:
    """

    files = []
    merge_conf_logs = []

    for occ in np.arange(self.params.validate.options.start_simul_occ,
                          self.params.validate.options.end_simul_occ,
                          self.params.validate.options.step_simulation):

        os.chdir(self.params.output.out_dir)

        # TODO Change filenames when Templating implemented #54
        #csv
        csv_filename = self.params.exhaustive.output + "_occ_{}_b_{}.csv".format(str(occ).replace(".","_"),
                                                                            str(self.params.validate.options.set_b).replace(
                                                                                ".","_"))
        files.append(csv_filename)
        assert os.path.exists(csv_filename), "CSV file {} missing".format(csv_filename)

        #png
        png_filename = self.params.exhaustive.output + "_occ_{}_b_{}.png".format(str(occ).replace(".","_"),
                                                                            str(self.params.validate.options.set_b).replace(
                                                                                ".", "_"))
        files.append(png_filename)
        assert os.path.exists(png_filename), "png file {} missing".format(png_filename)

        #simul pdb
        pdb_filename = self.params.input.xtal_name + "_refine_occ_{}.pdb".format(str(occ).replace(".","_"))
        files.append(pdb_filename)
        assert os.path.exists(pdb_filename), "pdb file {} missing".format(pdb_filename)

        #simul_pdb_set_b
        if self.params.validate.options.set_b is not None:

            set_b_pdb = pdb_filename.rstrip(".pdb") + "_set_b_{}".format(str(self.params.validate.options.set_b).replace(
                                                                                ".", "_"))
            files.append(set_b_pdb)
            assert os.path.exists(set_b_pdb), "pdb file {} missing".format(set_b_pdb)
            pass

            mtz_filename = set_b_pdb + ".mtz"
        else:
            mtz_filename = pdb_filename + ".mtz"

        #mtz
        files.append(mtz_filename)
        assert os.path.exists(mtz_filename), "mtz file {} missing".format(mtz_filename)

        if self.params.validate.options.use_qsub():
            pass

        # merge-conf logs
        merge_conf_log = "_occ_{}_b_{}_merge_conformations.log".format(str(occ).replace(".","_"),
                                                   str(self.params.validate.options.set_b).replace(".", "_"))
        merge_conf_logs.append(merge_conf_log)
        assert os.path.exists(merge_conf_log), "Merge conformation log {} missing".format(merge_conf_log)

    logs = [f for f in os.listdir(os.path.join(self.params.output.out_dir,self.params.output.log_dir))
            if isfile(os.path.join(self.params.output.out_dir,self.params.output.log_dir, f))]

    # es logs
    es_logs =[]
    for log in logs:
        if log.startswith(self.params.exhaustive.output.log_name):
            es_logs.append(log)

    assert len(es_logs)== len(np.arange(self.params.validate.options.start_simul_occ,
                          self.params.validate.options.end_simul_occ,
                          self.params.validate.options.step_simulation)), "Incorrect number of " \
                                                                     "Exhaustive Search Logs:\n{}\n".format(es_logs)

    # validate log
    validate_logs(log)
    for log in logs:
        if log.startswith(self.params.validate.output.log_name):
            validate_logs.append(log)

    assert(len(validate_log) == 1), "Incorrect number of validation logs: \n{}\n".format(validate_log)

    # summary png
    summary_png = self.params.input.xtal_name()

    assert os.path.exists("multi-state-restraints.refmac.self.params"), "Missing: multi-state-restraints.refmac.params"
    files.append("multi-state-restraints.refmac.params")
    assert os.path.exists("multi-state-restraints.phenix.params"), "Missing: multi-state-restraints.phenix.params"
    files.append("multi-state-restraints.refmac.params")

    # Check no other files are generated
    generated_files = [f for f in os.listdir(self.params.output.out_dir) if isfile(os.path.join(self.params.output.out_dir, f))]


    assert len(files) == len(generated_files), "There are {} extra files in the output folder".format(len(files) -
                                                                                                      len(generated_files))

if __name__ == '__main__':
    unittest.main()
