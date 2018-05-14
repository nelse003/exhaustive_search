import os

import numpy as np

input_data_mtz = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/refine.mtz"

for occupancy in np.arange(0,1.01,0.05):

    input_simul_mtz = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation.py/simul_{}.mtz".format(
        str(occupancy).replace(".", "_"))
    output_simul_mtz = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation.py/simul_cad_{}.mtz".format(
        str(occupancy).replace(".", "_"))

    cmd = ('cad hklin1 {} hklin2 {} hklout {} << eof\n'.format(input_simul_mtz, input_data_mtz,output_simul_mtz) +
           '  monitor BRIEF\n'
           '  labin file 1 E1=H E2=K E3=L E4=FOBS(+) E5=SIGFOBS(+) E6=FOBS(-) E7=SIGFOBS(-)\n'
           #'  labout file 1 E1=H E2=K E3=L E4=FOBS(+) E5=SIGFOBS(+) E6=FOBS(-) E7=SIGFOBS(-)\n'
           '  labin file 2 E1=PHIC \n'
           #'  labout file 2 E1=PHIC\n'
           'eof\n')

    print(cmd)

    os.system(cmd)

    exit()