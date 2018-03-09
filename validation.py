import os

import numpy as np


def buffer_validation():

    input_pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/refine.pdb"
    input_mtz = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/refine.mtz"
    out_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/buffer_vary/"
    xtal_name = "NUDT7A-x1237"

    for buffer in np.arange(0.75,30.01,0.25):

        out_dir = os.path.join(out_path, "{}_{}".format(xtal_name, buffer))
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        buffer_txt = str(buffer).replace(".","_")
        sh_file = "{}_{}.sh".format(xtal_name,buffer_txt)
        #
        # file = open(os.path.join(out_dir, sh_file),'w')
        #
        # file.write("#!/bin/bash\n")
        # file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
        # file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
        # file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive_search.py "\
        #            "input.pdb={} input.mtz={} options.buffer={} "\
        #            "xtal_name={} output.out_dir={} \n".format(input_pdb, input_mtz, buffer,xtal_name, out_dir))
        #
        # file.close()

        os.system("qsub -o {} -e {} {}".format(os.path.join(out_dir,"output.txt"),
                                             os.path.join(out_dir,"error.txt")
                                            ,os.path.join(out_dir, sh_file)))
        #scatter_plot(os.path.join(out_dir,"u_iso_occupancy_vary"))

