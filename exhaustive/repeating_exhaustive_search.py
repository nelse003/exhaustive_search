from __future__ import print_function

import os
import sqlite3
import sys

import libtbx.phil
import pandas as pd

from exhaustive.exhaustive import run as exhaustive_search
from ..phil import master_phil
##############################################################

PROGRAM = 'Repeat Exhaustive Search'
DESCRIPTION = """
    Take in csv of pdb and mtz locations,
    run a exhaustive search of B factor and occupancy to determine unique minima.
"""
blank_arg_prepend = {'.pdb': 'pdb=','.mtz': 'mtz=','.csv': 'csv='}

########################################################################
import logging
import datetime

logging.basicConfig(filename=datetime.datetime.now().strftime('/dls/science/groups/i04-1/elliot-dev/Work/' \
                                                              'exhaustive_search/logs/exhaustive_search_%H_%M_%d_%m_%Y.log'),
                    level=logging.DEBUG)
logging = logging.getLogger(__name__)
########################################################################



        #



    #
    #     #### For Plotting ####
    #
    #     output_folder = "test_runs/{}".format(xtal_name)
    #     output_path = os.path.join(os.getcwd(), output_folder)
    #     output_path_base = os.path.join(os.getcwd(), "output_DCP2_refinements")
    #
    #     if not os.path.exists(output_path_base):
    #         os.mkdir(output_path_base)
    #
    #     if not os.path.exists(output_path):
    #         os.mkdir(output_path)
    #     os.chdir(output_folder)
    #     csv_name = 'u_iso_occupancy_vary'
    #     print(os.getcwd())
    #     os.rename(csv_name, csv_name + ".csv")
    #     scatter_plot(csv_name)
            #     os.chdir("../../"


if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)


