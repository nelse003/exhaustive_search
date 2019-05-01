"""Exhaustive search.

Take in pdb and mtz, or csv of pdb and mtz locations,
run a exhaustive search of B factor and occupancy to
determine unique minima of mean(|Fo-Fc|).
"""

#  Imports
from __future__ import division
from __future__ import print_function

import logging
import os
import sys
import csv

from utils.phil import master_phil
from utils.utils_ccp4 import iter_u_iso_occ

from xtal_model_data import XtalModelData
from occ_b_loop import OccBLoopCaller

##############################################################
PROGRAM = "Exhaustive Search"
DESCRIPTION = """
    Take in pdb and mtz, or csv of pdb and mtz locations,
    run a exhaustive search of B factor and occupancy to
    determine unique minima.
"""
blank_arg_prepend = {".pdb": "pdb=", ".mtz": "mtz=", ".csv": "csv="}
##############################################################

# TODO Split, and move running function move to bin

logger = logging.getLogger(__name__)


def run(params):
    """ 
    Minima of |Fo-Fc| of superposed model across Occupancy and B factor

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
    params: libtbx.phil.scope_extract
        python object from phil file,
        edited with any additional parameters

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

    if not os.path.exists(os.path.join(params.output.out_dir, params.output.log_dir)):
        os.makedirs(os.path.join(params.output.out_dir, params.output.log_dir))

    logger.info("Running Exhaustive Search \n\n")

    modified_phil = master_phil.format(python_object=params)

    logger.info("Current Parameters")
    logger.info(master_phil.format(python_object=params).as_str())
    logger.info("Parameters Different from default")
    logger.info(master_phil.fetch_diff(source=modified_phil).as_str())

    logger.info("{}: running exhaustive search".format(str(params.input.xtal_name)))

    # create a list of occupancies and b factors to loop over
    u_iso_occ = iter_u_iso_occ(params)

    # Generate XtalModelData from pdb and mtz
    xtal_model_data = XtalModelData(params)

    # Define a instance of the loop
    occ_b_loop = OccBLoopCaller(xtal_model_data=xtal_model_data)

    # Use map() to loop over supplied list of occupancies and b factors
    sum_fofc_results = map(occ_b_loop, u_iso_occ)

    logger.info(
        "Loop finished.\n"
        "Writing bound occupancy, ground_occupancy, u_iso, "
        "mean |Fo-Fc| to CSV: {}".format(params.exhaustive.output.csv_name)
    )

    with open(params.exhaustive.output.csv_name, "w") as f1:
        writer = csv.writer(f1, delimiter=",", lineterminator="\n")
        writer.writerows(sum_fofc_results)
        sys.stdout.flush()


if __name__ == "__main__":
    from giant.jiffies import run_default

    run_default(
        run=run,
        master_phil=master_phil,
        blank_arg_prepend=blank_arg_prepend,
        args=sys.argv[1:],
    )
