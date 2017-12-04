import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def scatter_plot(csv_name):
    # Load data from CSV
    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    y = (8*np.pi**2)*y**2

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    ax.scatter(x,y,z)

    plt.xlabel("Occupancy")
    plt.ylabel("B_iso")
    ax.set_zlabel("FoFc")

    plt.show()

    plt.savefig(csv_name)
    plt.close()

def bounded_2d_scatter(atom_name,lower_bound,upper_bound):
    # Load data from CSV
    data = np.genfromtxt('{}.csv'.format(atom_name), delimiter=',', skip_header=0)

    reduced_data = data[np.where(np.logical_and(data[:,2] >= lower_bound,data[:,2] <= upper_bound))]

    x = reduced_data[:,0]
    y = reduced_data[:,1]

    y = (8 * np.pi ** 2) * y ** 2

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlabel("Occupancy")
    plt.ylabel("B_iso")
    plt.title("Fofc values between {} and {}".format(lower_bound,upper_bound))

    ax.scatter(x, y)
    plt.savefig("{}_reduced_".format(atom_name))
    plt.close()

def colourbar_2d_scatter(atom_name):
    # Load data from CSV
    data = np.genfromtxt('{}.csv'.format(atom_name), delimiter=',', skip_header=0)

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    y = (8 * np.pi ** 2) * y ** 2

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.xlabel("Occupancy")
    plt.ylabel("B_iso")

    sc = ax.scatter(x, y, c=z, marker =".")

    plt.colorbar(sc)
    plt.savefig("{}_colorbar".format(atom_name))
    plt.close()

###########################################################

# Per atom plots
"""
for i in range(71,79):
    atom_name = "HETATM_3{}".format((3 - len(str(i))) * '0' + str(i))
#    scatter_plot(atom_name)
#    bounded_2d_scatter(atom_name,lower_bound = -0.05, upper_bound = 0.05)
    colourbar_2d_scatter(atom_name)

#bounded_2d_scatter("HETATM_3077",-0.3,0.3)
#bounded_2d_scatter("HETATM_3071",-0.3,0.3)
#bounded_2d_scatter("HETATM_3072",-0.3,0.3)

"""

# Per residue plot
print os.getcwd()
os.chdir("../output")
scatter_plot("mean_point_near_lig_fofc")
#bounded_2d_scatter("LIG",-0.1,0.1)