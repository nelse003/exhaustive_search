import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from exhaustive.utils.utils import get_fofc_from_csv, get_minimum_fofc, round_step, b_to_u_iso, u_iso_to_b_fac
from mpl_toolkits.mplot3d import Axes3D

def scatter_plot(csv_name, three_dim_plot=True ,title_text=None ):

    """ Scatter plots of occupancy, U_iso and mean |fo_fc| from csv output """

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

    if title_text is not None:
        plt.title(title_text)

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

def connectpoints(x,y,x_1,y_1,p1):

    """ Draw lines between two lists of points"""

    x1, x2 = x[p1], x_1[p1]
    y1, y2 = y[p1], y_1[p1]
    plt.plot([x1,x2],[y1,y2],'k-')

def connectpoint(x,y,x_1,y_1):

    """ Draw line between two points """

    plt.plot([x,x_1],[y,y_1],'k-')

def connectpoints_3d(x,y,z,x_1,y_1,z_1,p1):

    """ Draw lines between two sets od points in 3d"""

    x1, x2 = x[p1], x_1[p1]
    y1, y2 = y[p1], y_1[p1]
    z1, z2 = z[p1], z_1[p1]
    plt.plot([x1,x2],[y1,y2],[z1,z2],'k-')


def plot_3d_fofc_occ(start_occ, end_occ, step, dataset_prefix, set_b):

    """ Plot the difference in occupancy & mean(|fo-fc|) at the simulated occupancy and the minima. """

    min_fofcs = []
    min_occs = []
    min_b_facs = []
    fofcs = []
    occs = []
    b_facs = []

    for lig_occupancy in np.arange(start_occ, end_occ + (step / 5), step):
        csv_name = "occ_{}_b_{}_u_iso".format(str(lig_occupancy).replace(".","_"),str(set_b).replace(".","_"))
        min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_name)
        fofc = get_fofc_from_csv(csv_name,lig_occupancy, round_step(b_to_u_iso(set_b)), step)
        fofcs.append(fofc)
        occs.append(lig_occupancy)
        b_facs.append(set_b)
        min_b_facs.append(u_iso_to_b_fac(min_u_iso))
        min_fofcs.append(fo_fc_at_min)
        min_occs.append(min_occ)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    min_plot, = ax.plot(min_occs, min_b_facs, min_fofcs,'k+')
    occ_plot, = ax.plot(occs, b_facs, fofcs, 'ro')

    for i in np.arange(0, len(occs)):
        connectpoints_3d(occs,b_facs,fofcs,min_occs,min_b_facs,min_fofcs,i)

    ax.legend((min_plot,occ_plot),
              ('Minima of mean |Fo-Fc|','Mean |Fo-Fc| at simulated occupancy'),
              prop={"size": 8},
              numpoints=1,
              bbox_to_anchor=(1, 1),
              bbox_transform=plt.gcf().transFigure)

    ax.set_xlabel("Occupancy")
    ax.set_ylabel("B factor")
    ax.set_zlabel("Mean |Fo-Fc|")

    plt.title("{}: Delta mean|Fo-Fc| and Delta Occupancy".format(dataset_prefix), fontsize=10)
    plt.savefig("{}-3d-delta_fofc_occ.png".format(dataset_prefix))

def plot_fofc_occ(start_occ, end_occ, step, dataset_prefix, set_b):

    """ Plot the difference in occupancy/ fofc at the simulated occupancy and minima.  """

    min_fofcs = []
    min_occs = []
    fofcs = []
    occs = []

    for lig_occupancy in np.arange(start_occ, end_occ + (step / 5), step):
        csv_name = "occ_{}_b_{}_u_iso".format(str(lig_occupancy).replace(".","_"),str(set_b).replace(".","_"))
        min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_name)
        fofc = get_fofc_from_csv(csv_name,lig_occupancy, round_step(b_to_u_iso(set_b)), step)
        fofcs.append(fofc)
        occs.append(lig_occupancy)
        min_fofcs.append(fo_fc_at_min)
        min_occs.append(min_occ)

    fig, ax = plt.subplots()
    min_plot, = ax.plot(min_occs, min_fofcs,'k+')
    occ_plot, = ax.plot(occs, fofcs, 'ro')

    for i in np.arange(0, len(occs)):
        connectpoints(occs,fofcs,min_occs,min_fofcs,i)

    ax.legend((min_plot,occ_plot),
              ('Minima of mean |Fo-Fc|','Mean |Fo-Fc| at simulated occupancy'),
              prop={"size": 8},
              numpoints=1,
              bbox_to_anchor=(1, 1),
              bbox_transform=plt.gcf().transFigure)

    ax.set_xlabel("Occupancy")
    ax.set_ylabel("Mean |Fo-Fc|")

    plt.title("{}: Delta mean|Fo-Fc| and Delta Occupancy".format(dataset_prefix), fontsize=10)
    plt.savefig("{}-delta_fofc_occ.png".format(dataset_prefix))


