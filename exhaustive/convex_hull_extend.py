from exhaustive import run as exhaustive
from phil import master_phil

params =  master_phil.extract()
params.testing.testing = True
params.input.xtal_name = "test"

exhaustive(params=params)