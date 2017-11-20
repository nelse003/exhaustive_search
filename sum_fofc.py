import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#######################################################
def scatter_plot_supply_z(atom_name,z):
    # Load data from CSV
    data = np.genfromtxt('{}.csv'.format(atom_name), delimiter=',', skip_header=0)

    x = data[:,0]
    y = data[:,1]

    y = (8*np.pi**2)*y**2

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    ax.scatter(x,y,z)

    plt.xlabel("Occupancy")
    plt.ylabel("B_iso")
    ax.set_zlabel("mean of FoFc over ligand")

    plt.savefig("mean_fofc")
    plt.close()

def bounded_2d_scatter_supply_z(atom_name,lower_bound,upper_bound,z):
    # Load data from CSV
    data = np.genfromtxt('{}.csv'.format(atom_name), delimiter=',', skip_header=0)

    reduced_data = data[np.where(np.logical_and(z >= lower_bound,z <= upper_bound))]

    x = reduced_data[:,0]
    y = reduced_data[:,1]

    y = (8 * np.pi ** 2) * y ** 2

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlabel("Occupancy")
    plt.ylabel("B_iso")
    plt.title("Fofc values between {} and {}".format(lower_bound,upper_bound))

    ax.scatter(x, y)
    plt.savefig("mean_fofc_reduced_1".format(atom_name))
    plt.close()

#########################################################

os.chdir('/home/enelson/Dropbox/DPhil/exhaustive_search/output_a')
# Loop over all atoms in an altloc, sum mFo-DFc for all B_iso and Occ

# This is a very stupid way to select altloc C of the LIG
# TODO make this more general to select ligand
all_fofc=[]

for i in range(72,141,2):
    atom_name = "HETATM_3{}".format((3 - len(str(i))) * '0' + str(i))
    data = np.genfromtxt('{}.csv'.format(atom_name), delimiter=',', skip_header=0)
    fofc_values = data[:,2]
    all_fofc.append(fofc_values)

mean_fofc = sum(all_fofc)/len(range(72,141,2))

scatter_plot_supply_z("HETATM_3140",mean_fofc)
bounded_2d_scatter_supply_z("HETATM_3140",-0.1,0.1,mean_fofc)





