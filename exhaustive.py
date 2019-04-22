"""Exhaustive search.

Take in pdb and mtz, or csv of pdb and mtz locations,
run a exhaustive search of B factor and occupancy to
determine unique minima of mean(|Fo-Fc|).
"""

#  Imports
from __future__ import division
from __future__ import print_function

import datetime
import logging
import os
import sys

from utils.phil import master_phil

from xtal_model_data import XtalModelData

##############################################################
PROGRAM = 'Exhaustive Search'
DESCRIPTION = """
    Take in pdb and mtz, or csv of pdb and mtz locations,
    run a exhaustive search of B factor and occupancy to
    determine unique minima.
"""
blank_arg_prepend = {'.pdb': 'pdb=', '.mtz': 'mtz=', '.csv': 'csv='}
##############################################################

logger = logging.getLogger(__name__)

def run(params):
    """ 

    Load in pdb and mtz file.
    Process pdb and mtz to produce:
     fmodel,
     iotbx protein hierarchy,
     xrs: xray structure
     inputs, a object used my mmtbx utils
     crystal gridding

    These are used to determine the minima of mean(|Fo-Fc|) over grid
    points defined by a search of Occupancy and B factor for atoms that
    are part of a bound ligand, or the changed protein that the bound
    atom is nearby.

    Parameters
    ----------------
    params: 
        A extracted python object from the master phil file. 
        This defines the settings of the settings and I/O. 

    Returns
    ----------------
    None:
        Function directly returns nothing
        
    Notes
    ----------------
    A CSV file containing the exhaustive search output:
    > Occupancy of Bound State
    > Occupancy of Ground State
    > U_iso of state (Isotropic B factor)
    > Mean |Fo-Fc| value over selected points
    
    is output. It is stored under filename provided in:
    params.exhaustive.options.csv_name
    
    """

    if not os.path.exists(params.output.out_dir):
        os.makedirs(params.output.out_dir)

    if not os.path.exists(os.path.join(params.output.out_dir,
                                       params.output.log_dir)):
        os.makedirs(os.path.join(
            params.output.out_dir,
            params.output.log_dir))

    log_time = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M.log")
    log_path = os.path.join(params.output.out_dir,
                            params.output.log_dir,
                            params.exhaustive.output.log_name + log_time)

    logger.info("Running Exhaustive Search \n\n")

    modified_phil = master_phil.format(python_object=params)

    logger.info("Current Parameters")
    logger.info(master_phil.format(python_object=params).as_str())
    logger.info("Parameters Different from default")
    logger.info(master_phil.fetch_diff(source=modified_phil).as_str())

    logger.info("{}: running exhaustive search".format(str(params.input.xtal_name)))

    xtal_model_data = XtalModelData(params)

    logger.info("Organising output directory")
    os.chdir(params.output.out_dir)

    logger.info("Run main calculation of |Fo-Fc| at grid points near ligand")

    try:
        xtal_model_data.calculate_mean_fofc()

    except UnboundLocalError:
        raise

    os.chdir("../../")


if (__name__ == "__main__"):
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil,
                blank_arg_prepend=blank_arg_prepend, args=sys.argv[1:])


