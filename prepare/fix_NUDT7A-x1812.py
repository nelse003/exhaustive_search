import os

pandda_input = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1812/dimple.pdb"
pandda_model = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1812/NUDT7A-x1812-pandda-model.pdb"

# os.mkdir("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/NUDT7A-x1812-new-test")
os.chdir("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/NUDT7A-x1812-new-test")

os.system("giant.merge_conformations major={} minor={}".format(pandda_input, pandda_model))
