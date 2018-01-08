import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def scatter_plot(csv_name, three_dim_plot=True):
    # Load data from CSV
    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)

    if len(data[0]) == 3:
        occ = data[:,0]
        u_iso = data[:,1]
        fo_fc = data[:,2]

    if len(data[0]) == 4:
        occ = data[:, 0]
        u_iso = data[:, 2]
        fo_fc = data[:, 3]

    b_iso = (8*np.pi**2)*u_iso**2

    fig = plt.figure()

    if three_dim_plot:
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(occ, b_iso, fo_fc)
        plt.xlabel("Occupancy")
        plt.ylabel("B_iso")
        ax.set_zlabel("Fo-Fc")
    else:
        ax = fig.add_subplot(111)
        ax.scatter(occ, fo_fc)

        plt.xlabel("Occupancy")
        plt.ylabel("Fo-Fc")

    plt.savefig(csv_name)
    plt.close()


def scatter_plot_8col(csv_name):
    # Load data from CSV
    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)

    lig_occ = data[:,0]
    bound_occ = data[:,1]
    u_iso = data[:,2]
    mean_local_abs_fofc = data[:,3]
    mean_abs_fofc_value = data[:,4]
    mean_fofc = data[:,5]
    mean_grid_fofc = data[:,6]
    mean_grid_abs_fofc = data[:,7]


    B_iso = (8*np.pi**2)*u_iso**2

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    ax.scatter(lig_occ, B_iso, mean_local_abs_fofc - mean_grid_fofc)

    plt.xlabel("Lig Occupancy")
    plt.ylabel("B iso")
    ax.set_zlabel("mean Fo - Fc")

    #plt.show()

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
# TODO make work
def get_steps(csv_name):
    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)

    occ = data[:, 0]
    u_iso = data[:, 2]
    fo_fc = data[:, 3]

    print occ

    smallest_occ_step =[x-occ[i-1] for i, x in enumerate(occ)][1:]
    print smallest_occ_step


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

#os.chdir("../output_0A_frac_fix")

# Per residue plot
#print os.getcwd()

#scatter_plot_4col("mean_ground_bound_5A_buffer")
#scatter_plot_2d("fixed_B_5A_covary")
#print os.getcwd()
#os.chdir("../output_4kjt")
#scatter_plot_4col("0A_bound_ground_covary_frac_fix")
#bounded_2d_scatter("LIG",-0.1,0.1)

