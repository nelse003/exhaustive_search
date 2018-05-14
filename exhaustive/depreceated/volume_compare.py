import numpy as np
import scipy.spatial as ss

from exhaustive.process.minima import get_minimum_fofc


def trapezoidal_area(xyz):
    """Calculate volume under a surface defined by irregularly spaced points
    using delaunay triangulation. "x,y,z" is a <numpoints x 3> shaped ndarray."""
    d = ss.Delaunay(xyz[:,:2])
    tri = xyz[d.vertices]

    a = tri[:,0,:2] - tri[:,1,:2]
    b = tri[:,0,:2] - tri[:,2,:2]
    proj_area = np.cross(a, b).sum(axis=-1)
    zavg = tri[:,:,2].sum(axis=1)
    vol = zavg * np.abs(proj_area) / 6.0
    return vol.sum()

start_occ = 0.05
end_occ = 0.95
step = 0.05

all_fofc_gaps = []
for simul_occ in np.arange(start_occ,end_occ+step/5,step):

    csv_path = ("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation.py/validation_bound_ground/" +
                "occ_{}_b_40_u_iso".format(str(simul_occ).replace('.','_')))

    data = np.genfromtxt('{}.csv'.format(csv_path), delimiter=',', skip_header=0)
    _,_,min_fofc = get_minimum_fofc(csv_path)
    fofc_gaps = data[:,3] - min_fofc
    sum_fofc_gaps = np.sum(fofc_gaps)
    all_fofc_gaps.append(sum_fofc_gaps)

mean_fofc_box = (np.mean(all_fofc_gaps)*0.15)/len(data)
print(mean_fofc_box)
print("u_iso:{}".format(np.sqrt(40/(8 * np.pi ** 2))))

for simul_occ in np.arange(start_occ,end_occ+step/5,step):

    csv_path = ("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation.py/validation_bound_ground/" +
                "occ_{}_b_40_u_iso".format(str(simul_occ).replace('.','_')))

    data = np.genfromtxt('{}.csv'.format(csv_path), delimiter=',', skip_header=0)
    _,_,min_fofc = get_minimum_fofc(csv_path)

    fofc_box = min_fofc + mean_fofc_box

    print("Simulated Occupancy {}".format(simul_occ))
    print(len(data[data[:,3] < fofc_box]))

    print(np.sort(data[data[:,3] < fofc_box]))

    print("________________________________________________")

