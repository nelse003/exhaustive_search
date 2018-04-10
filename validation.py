import os

import matplotlib.pyplot as plt
import numpy as np
from iotbx.pdb import hierarchy

from plot_exhaustive_search import scatter_plot
from process_exhaustive_search import get_minimum_fofc
from select_occupancy_groups import get_occupancy_groups


def buffer_validation():

    input_pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/refine.pdb"
    input_mtz = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/refine.mtz"
    out_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/buffer_vary/"
    xtal_name = "NUDT7A-x1237"

    for buffer in np.arange(0.0,0.76,0.25):

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

        # os.system("qsub -o {} -e {} {}".format(os.path.join(out_dir,"output.txt"),
        #                                      os.path.join(out_dir,"error.txt")
        #                                     ,os.path.join(out_dir, sh_file)))
        scatter_plot(os.path.join(out_dir,"u_iso_occupancy_vary"))


def occ_loop_simulate_exp_data():

    input_pdb = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/refine.pdb"
    input_mtz = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1237/NUDT7A-x1237.mtz"
    xtal_name = "NUDT7A-x1237"

    os.chdir("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation")

    for occupancy in np.arange(0,1.01,0.05):
        # pdb_in = hierarchy.input(input_pdb)
        #
        # for chain in pdb_in.hierarchy.only_model().chains() :
        #     for residue_group in chain.residue_groups() :
        #         for atom_group in residue_group.atom_groups() :
        #             if atom_group.resname == "LIG":
        #                 for atom in atom_group.atoms() :
        #                     atom.set_occ(occupancy)
        #
        occ_pdb = "NUDT7A-x1237-occ_{}.pdb".format(str(occupancy).replace(".","_"))
        # f = open(occ_pdb, "w")
        # f.write(pdb_in.hierarchy.as_pdb_string(
        #     crystal_symmetry=pdb_in.input.crystal_symmetry()))
        # f.close()
        #
        log_file = "simul_from_refine_log_{}".format(str(occupancy).replace(".","_"))
        mtz_out = "simul_occ_{}.mtz".format(str(occupancy).replace(".","_"))

        os.system("ccp4-python ../simulate_experimental_data.py input.xray_data.file_name={} "
                  "model.file_name={} input.xray_data.label=\"I(+),SIGI(+),I(-),SIGI(-),merged\" "
                  "output.logfile={} output.hklout={}".format(input_mtz, occ_pdb, log_file,mtz_out))

        # input_mtz = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/simul_{}.mtz".format(str(occupancy).replace(".","_"))
        #
        # out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/"
        # csv_name = "occ_u_iso_{}".format(str(occupancy).replace(".","_"))
        #
        # sh_file = "{}_{}.sh".format(xtal_name, occupancy)
        #
        # file = open(os.path.join(out_dir, sh_file),'w')
        #
        # file.write("#!/bin/bash\n")
        # file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
        # file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
        #
        # file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive_search.py"
        #            " input.pdb={} input.mtz={} output.out_dir={} xtal_name={} "
        #            "options.csv_name={}".format(input_pdb,input_mtz,out_dir,xtal_name,csv_name))
        # file.close()
        #
        # os.system("qsub -o {} -e {} {}".format(os.path.join(out_dir,"output_{}.txt".format(str(occupancy).replace(".","_"))),
        #                                        os.path.join(out_dir,"error_{}.txt".format(str(occupancy).replace(".","_"))),
        #                                        os.path.join(out_dir, sh_file)))

    # args = [input_pdb, input_mtz]
    #
    # pdb_inp = iotbx.pdb.input(input_pdb)
    #
    # inputs = mmtbx.utils.process_command_line_args(args = args)
    # ph = pdb_inp.construct_hierarchy()
    # xrs = ph.extract_xray_structure(
    #     crystal_symmetry = inputs.crystal_symmetry)
    # reflection_files = inputs.reflection_files
    # rfs = reflection_file_utils.reflection_file_server(
    #     crystal_symmetry=inputs.crystal_symmetry,
    #     force_symmetry=True,
    #     reflection_files=reflection_files,
    #     err=StringIO())
    # determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
    #     reflection_file_server=rfs,
    #     keep_going=True,
    #     log=StringIO())
    # r_free_flags = determined_data_and_flags.r_free_flags
    # f_obs = determined_data_and_flags.f_obs
    # crystal_gridding = f_obs.crystal_gridding(
    #     d_min             = f_obs.d_min(),
    #     symmetry_flags    = maptbx.use_space_group_symmetry,
    #     resolution_factor = 1./4)
    #
    # mask_params = mmtbx.masks.mask_master_params.extract()
    # mask_params.ignore_hydrogens=False
    # mask_params.ignore_zero_occupancy_atoms=False
    # fmodel = mmtbx.f_model.manager(
    #     r_free_flags   = r_free_flags
    #     mask_params    = mask_params,
    #     xray_structure = xrs)
    # fmodel.update_all_scales()
    #
    # fcfc_map, fcfc  = compute_maps(fmodel=fmodel,crystal_gridding=crystal_gridding,map_type="mFo-DFc")
    # mtz_dataset = fcfc.as_mtz_dataset(column_root_label="FOFCWT")
    # mtz_object = mtz_dataset.mtz_object()
    # mtz_object.write(file_name="testing_{}_{}.mtz".format(bound_occupancy, u_iso))

def refmac_0_cyc(input_mtz, input_pdb, output_pdb, output_mtz, occupancy):

    with open("refmac_0_cyc_occ_{}.sh".format(str(occupancy).replace(".", "_")), 'w') as f:
        f.write("#!/bin/bash \n")
        f.write(
            "source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
        f.write("refmac5 HKLIN {} \\\n".format(input_mtz))
        f.write("HKLOUT {} \\\n".format(output_mtz))
        f.write("XYZIN {} \\\n".format(input_pdb))
        f.write("XYZOUT {} \\\n".format(output_pdb))
        f.write("LIBIN /dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model/NUDT7A-x1740/NUOOA000181a.cif \\\n")
        f.write("LIBOUT /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/NUOOA000181a.cif \\\n")
        f.write(" << EOF > refmac_{}.log".format(str(occupancy).replace(".", "_")))
        f.write("""
make -
    hydrogen ALL -
    hout NO -
    peptide NO -
    cispeptide YES -
    ssbridge YES -
    symmetry YES -
    sugar YES -
    connectivity NO -
    link YES
refi -
    type REST -
    resi MLKF -
    meth CGMAT -
    bref ISOT
ncyc 0
scal -
    type SIMP -
    LSSC -
    ANISO -
    EXPE
weight AUTO
solvent YES
monitor MEDIUM -
    torsion 10.0 -
    distance 10.0 -
    angle 10.0 -
    plane 10.0 -
    chiral 10.0 -
    bfactor 10.0 -
    bsphere 10.0 -
    rbond 10.0 -
    ncsr 10.0
labin  FP=FOBS SIGFP=SIGFOBS
labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT FOM=FOM
DNAME NUDT7A-x1851
END
EOF
    """)

    os.system("qsub refmac_0_cyc_occ_{}.sh".format(str(occupancy).replace(".", "_")))

def occ_loop_merge_confs_simulate(bound_state_pdb_path,
                                             ground_state_pdb_path,
                                             input_mtz,
                                             dataset_prefix,
                                             out_path,
                                             set_b = None):

    for lig_occupancy in np.arange(0.05, 0.96, 0.05):
        merged_pdb = os.path.join(out_path,
                                  "{}_refine_occ_{}.pdb".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        os.system("giant.merge_conformations input.major={} input.minor={} "
                  "major_occupancy={} minor_occupancy={} output.pdb={}".format(
            ground_state_pdb_path, bound_state_pdb_path, str(1 - lig_occupancy), str(lig_occupancy), merged_pdb))

        simulate_log = os.path.join(out_path,"{}_simul_{}.log".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))
        simulate_mtz = os.path.join(out_path,"{}_simul_{}.mtz".format(dataset_prefix, str(lig_occupancy).replace(".", "_")))

        if set_b is not None:
            merged_file_name, _ = os.path.splitext(merged_pdb)

            set_u_iso_all_occupancy_groups(input_pdb = merged_pdb,
                                           output_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_")),
                                           b_fac = set_b)

            merged_pdb = merged_file_name + "_set_b_{}.pdb".format(
                                               str(set_b).replace(".", "_"))

        os.system("ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/simulate_experimental_data.py "
                  "input.xray_data.file_name={} "
                  "model.file_name={} input.xray_data.label=\"F,SIGF\" "
                  "output.logfile={} output.hklout={}".format(input_mtz, merged_pdb, simulate_log, simulate_mtz))

        # Exhaustive search
        sh_file = "{}_occ_{}_b_{}.sh".format(dataset_prefix,
                                             str(lig_occupancy).replace(".", "_"),
                                             str(set_b).replace(".", "_"))
        csv_name = "occ_{}_b_{}_u_iso".format(str(lig_occupancy).replace(".", "_"),str(set_b).replace(".", "_"))

        with open(os.path.join(out_path, sh_file),'w') as file:

            file.write("#!/bin/bash\n")
            file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
            file.write("source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")

            file.write("$CCP4/bin/ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive_search.py"
                       " input.pdb={} input.mtz={} output.out_dir={} xtal_name={} "
                       "options.csv_name={}".format(merged_pdb,simulate_mtz,out_path,dataset_prefix,csv_name))

        os.system("qsub -o {} -e {} {}".format(os.path.join(out_path,"output_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                               os.path.join(out_path,"error_{}.txt".format(str(lig_occupancy).replace(".","_"))),
                                               os.path.join(out_path, sh_file)))

        #Refmac 0 cycles
        output_pdb = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.pdb".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        output_mtz = os.path.join(out_path,"{}_occ_{}_b_{}_refmac_0_cyc.mtz".format(dataset_prefix,
                                                                           str(lig_occupancy).replace(".","_"),
                                                                           str(set_b).replace(".", "_")))
        refmac_0_cyc(input_mtz = simulate_mtz, input_pdb = merged_pdb,
                     output_pdb = output_pdb , output_mtz = output_mtz,
                     occupancy= lig_occupancy)

def set_u_iso_all_occupancy_groups(input_pdb, output_pdb, b_fac):

    pdb_inp = hierarchy.input(input_pdb)
    occ_group = get_occupancy_groups(pdb= input_pdb)
    for group in occ_group[0]:
        for residue in group:

            #TODO Replace many for loops with better structure

            for chain in pdb_inp.hierarchy.only_model().chains():
                if chain.id == residue.get('chain'):
                    for residue_group in chain.residue_groups():
                        if residue_group.resseq == residue.get('resseq'):
                            for atom_group in residue_group.atom_groups():
                                for atom in atom_group.atoms() :
                                    atom.set_b(b_fac)
                                    print("Changed: {}".format(residue))
    with open(output_pdb,"w") as f:
        f.write(pdb_inp.hierarchy.as_pdb_string(crystal_symmetry=pdb_inp.input.crystal_symmetry()))

def read_occ_csv(csv_name):

    """Return numpy arrays from  a csv with lig_occ, ground_occ, u_iso, mean(|Fo-Fc|)"""

    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)
    occ = data[:, 0]
    u_iso = data[:, 1]
    fo_fc = data[:, 2]

    return occ, u_iso, fo_fc

def get_delta_occ(csv_name, occupancy):

    """ Find the difference between minima occupancy and supplied occupancy.
     
     Given a csv with lig_occ, ground_occ, u_iso, mean(|Fo-Fc|) csv and occupancy value.
     """

    min_occ, _, _ = get_minimum_fofc(csv_name)
    return occupancy - min_occ

def get_delta_u_iso(csv_name):

    """ Find the minima u_iso and supplied u_iso.

     Given a csv with lig_occ, ground_occ, u_iso, mean(|Fo-Fc|) csv and u_iso value. 
     Useful when using a set B factor.
     """

    _, min_u_iso, _ = get_minimum_fofc(csv_name)
    return u_iso - min_u_iso

def b_to_u_iso(b_fac):

    u_iso = np.sqrt(b_fac/(8 * np.pi ** 2))
    return u_iso

def u_iso_to_b_fac(u_iso):

    b_iso = (8 * np.pi ** 2) * u_iso ** 2
    return b_iso

def round_step(x, prec=2, base=.05):
  return round(base * round(float(x)/base),prec)

def get_delta_fofc(csv_name, occupancy, u_iso, step=0.05):

    """ Find the difference between mean(|Fo-Fc|) at minima and at supplied u_iso and occ"""

    fo_fc = get_fofc_from_csv(csv_name,occupancy, u_iso)

    min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_name)
    delta_fofc = fo_fc - fo_fc_at_min

    return delta_fofc

def get_fofc_from_csv(csv_name,occupancy, u_iso, step=0.05):

    occupancy = round_step(occupancy,base=step)
    u_iso = round_step(u_iso,base=step)

    data = np.genfromtxt('{}.csv'.format(csv_name), delimiter=',', skip_header=0)
    e=0.0001
    data_line = data[((occupancy-e)<data[:,0]) & ((occupancy+e)>data[:,0]) & ((u_iso-e)<data[:,2]) & ((u_iso+e)>data[:,2]) ]
    fo_fc = data_line[0][3]

    return fo_fc

def get_delta_fofc_over_occupancies(set_b,start_occ, end_occ, step = 0.05):

    """ Get delta occupancy across occupancies"""

    occupancies = []
    delta_fofcs = []

    for lig_occupancy in np.arange(start_occ, end_occ+(step/5), step):
        delta_fofc = get_delta_fofc("occ_{}_b_{}_u_iso".format(
                                    str(lig_occupancy).replace(".","_"),
                                    str(set_b).replace(".","_")),
                                    lig_occupancy,
                                    b_to_u_iso(set_b),
                                    step = step)
        occupancies.append(lig_occupancy)
        delta_fofcs.append(delta_fofc)

    return np.column_stack((occupancies,delta_fofcs))

def connectpoints(x,y,x_1,y_1,p1):
    x1, x2 = x[p1], x_1[p1]
    y1, y2 = y[p1], y_1[p1]
    plt.plot([x1,x2],[y1,y2],'k-')

def occ_diff(occs, min_occs):

    eps = 1e-15
    occ_diff = np.abs(np.array(occs) - np.array(min_occs))
    occ_diff[np.abs(occ_diff)< eps] = 0
    frac_zeroes = float(len(occ_diff) - np.count_nonzero(occ_diff))/float(len(occ_diff))
    mean_occ_diff = np.mean(occ_diff)

    return frac_zeroes,mean_occ_diff

def plot_fofc_occ(start_occ, end_occ, step, set_b):

    min_fofcs = []
    min_occs = []
    fofcs = []
    occs = []

    for lig_occupancy in np.arange(start_occ, end_occ + (step / 5), step):
        csv_name = "occ_{}_b_{}_u_iso".format(str(lig_occupancy).replace(".","_"),str(set_b).replace(".","_"))
        min_occ, min_u_iso, fo_fc_at_min = get_minimum_fofc(csv_name)
        fofc = get_fofc_from_csv(csv_name,lig_occupancy, round_step(b_to_u_iso(40)), step)
        fofcs.append(fofc)
        occs.append(lig_occupancy)
        min_fofcs.append(fo_fc_at_min)
        min_occs.append(min_occ)

    frac_zeroes, mean_occ_diff = occ_diff(occs, min_occs)

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

    plt.title("Delta mean|Fo-Fc| and Delta Occupancy", fontsize=10)
    plt.savefig("NUDT7A-x1740-set-b-delta_fofc_occ.png")


in_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/validation_bound_ground"
bound_state_pdb_path = os.path.join(in_path,"refine.output.bound-state.pdb")
ground_state_pdb_path = os.path.join(in_path,"refine.output.ground-state.pdb")
input_mtz = os.path.join(in_path,"NUDT7A-x1740.free.mtz")
dataset_prefix = "NUDT7A-x1740"
out_path = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/validation/validation_bound_ground"

#delta_fofc_array = get_delta_fofc_over_occupancies(40,0.05, 0.95, step = 0.05)
#print(delta_fofc_array)

plot_fofc_occ(0.05, 0.95, step = 0.05, set_b = 40)

# occ_loop_merge_confs_simulate(bound_state_pdb_path,
#                                          ground_state_pdb_path,
#                                          input_mtz,
#                                          dataset_prefix,
#                                          out_path,
#                                          set_b = 40)