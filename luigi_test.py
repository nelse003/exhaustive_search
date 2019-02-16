import luigi
import os
from parse_xchemdb import main as parse_xchemdb
from occ_xchemdb_histogram import main as plot_occ
class ParseXChemDB(luigi.Task):
    def requires(self):
        return None
    def output(self):
        return luigi.LocalTarget('/dls/science/groups/i04-1/'
                                 'elliot-dev/Work/'
                                 'exhaustive_parse_xchem_db/'
                                 'log_pdb_mtz.csv')
    def run(self):
        parse_xchemdb()

class OccConvergence(luigi.Task):
    def requires(self):
        return ParseXChemDB()
    def output(self):
        return luigi.LocalTarget("/dls/science/groups/i04-1/elliot-dev/"
                                 "Work/exhaustive_parse_xchem_db/"
                                 "occ_conv.csv")
    def run(self):
        os.system("ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/"
                  "exhaustive_search/occ_convergence_on_db.py")

class PlottingOccHistogram(luigi.Task):
    def requires(self):
        return OccConvergence()
    def output(self):
        return luigi.LocalTarget("/dls/science/groups/i04-1/elliot-dev/"
                                 "Work/exhaustive_parse_xchem_db/"
                                 "bound_occ_hist.png")
    def run(self):
        plot_occ()



if __name__ == '__main__':
    luigi.build([PlottingOccHistogram()], local_scheduler=True, workers=10)