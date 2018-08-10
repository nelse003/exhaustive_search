from exhaustive import run as exhaustive
from phil import master_phil

params =  master_phil.extract()
params.testing = True
exhaustive(params=params)