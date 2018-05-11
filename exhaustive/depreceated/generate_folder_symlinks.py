from __future__ import print_function
import os

input_initial_model_folder = "/dls/labxchem/data/2017/lb18145-3/processing/analysis/initial_model"
output_initial_model_folder = "/dls/labxchem/data/2017/lb18145-3/processing/Elliot/processing/analysis/initial_model"

for folder in os.listdir(input_initial_model_folder):
    if not os.path.exists(os.path.join(output_initial_model_folder,folder)):
        os.mkdir(os.path.join(output_initial_model_folder,folder))

    if os.path.exists(os.path.join(input_initial_model_folder,folder,"dimple.pdb")):
        os.symlink(os.path.join(input_initial_model_folder,folder,"dimple.pdb"),
                   os.path.join(output_initial_model_folder,folder,"dimple.pdb"))

    if os.path.exists(os.path.join(input_initial_model_folder,folder,"dimple.mtz")):
        os.symlink(os.path.join(input_initial_model_folder,folder,"dimple.mtz"),
                   os.path.join(output_initial_model_folder,folder,"dimple.mtz"))