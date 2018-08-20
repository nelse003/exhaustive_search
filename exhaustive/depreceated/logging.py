import logging
import datetime
from phil import master_phil

def start_exhaustive_logging(params):
    """Prepare logging.

    NOT CURRENTLY USED

    Logging for exhaustive search using python logging module.
    Includes a description of parameters, and parameter different to default.
    """
    log_time = datetime.datetime.now().strftime("_%Y_%m_%d_%H_%M.log")
    log_path = os.path.join(params.output.out_dir,
                            params.output.log_dir,
                            params.exhaustive.output.log_name + log_time)
    hdlr = logging.FileHandler(log_path)
    logging = logging.getLogger(__name__)
    formatter = logging.Formatter('%(asctime)s %(levelname)s \n %(message)s')
    hdlr.setFormatter(formatter)
    logging.addHandler(hdlr)
    logging.setLevel(0)
    logging.info("Running Exhaustive Search \n\n")

    modified_phil = master_phil.format(python_object=params)
    logging.info("Current Parameters")
    logging.info(master_phil.format(python_object=params).as_str())
    logging.info("Parameters Different from default")
    logging.info(master_phil.fetch_diff(source=modified_phil).as_str())

    return logging