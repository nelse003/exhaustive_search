import os

def refmac_0_cyc(input_mtz, input_pdb, output_pdb, output_mtz, input_cif, output_cif, occupancy):

    with open("refmac_0_cyc_occ_{}.sh".format(str(occupancy).replace(".", "_")), 'w') as f:
        f.write("#!/bin/bash \n")
        f.write(
            "source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
        f.write("refmac5 HKLIN {} \\\n".format(input_mtz))
        f.write("HKLOUT {} \\\n".format(output_mtz))
        f.write("XYZIN {} \\\n".format(input_pdb))
        f.write("XYZOUT {} \\\n".format(output_pdb))
        f.write("LIBIN {} \\\n".format(input_cif))
        f.write("LIBOUT {} \\\n".format(output_cif))
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

def submit_exhasutive_with_refmac_0(dataset_prefix, out_path, set_b = None):

    """ Qsub submission of sh files from occ_loop_merge_confs_simulate_with_refmac_0"""

    for lig_occupancy in np.arange(0.05, 0.96, 0.05):

        sh_file = "{}_occ_{}_b_{}.sh".format(dataset_prefix,
                                             str(lig_occupancy).replace(".", "_"),
                                             str(set_b).replace(".", "_"))

        os.system("qsub -o {} -e {} {}".format(
            os.path.join(out_path, "output_{}.txt".format(str(lig_occupancy).replace(".", "_"))),
            os.path.join(out_path, "error_{}.txt".format(str(lig_occupancy).replace(".", "_"))),
            os.path.join(out_path, sh_file)))

def refine_after_exhasutive_search(input_pdb, input_mtz, input_cif, refine_params, dataset_prefix, working_dir = None):

    os.chdir(working_dir)
    sh_file = "{}_quick_refine_exhaustive_search_minima.sh".format(dataset_prefix)
    out_prefix = "{}_refine_after_exhaustive_search".format(dataset_prefix)
    dir_prefix = out_prefix
    quick_refine_qsub(input_pdb=input_pdb, input_mtz=input_mtz, input_cif=input_cif, refine_params=refine_params,
                      sh_file=sh_file, out_prefix=out_prefix, dir_prefix=dir_prefix)


def quick_refine_qsub(input_pdb, input_mtz, input_cif, refine_params,
                      sh_file, out_prefix=None, dir_prefix=None):

    with open(os.path.join(out_path, sh_file), 'w') as file:
        file.write("#!/bin/bash\n")
        file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
        file.write(
            "source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
        file.write("giant.quick_refine input.pdb={} input.mtz={} input.cif={} "
                   "input.params={} output.out_prefix={} output.dir_prefix={} ".format(input_pdb, input_mtz,
                                                                                       input_cif, refine_params,
                                                                                       out_prefix, dir_prefix))

    os.system("qsub -o {} -e {} {}".format(
            os.path.join(out_path, "output_{}.txt".format(str(out_prefix))),
            os.path.join(out_path, "error_{}.txt".format(str(out_prefix))),
            os.path.join(out_path, sh_file)))


def quick_refine_repeats(start_occ, end_occ,step, dataset_prefix, set_b, out_path, input_cif):

    params = "multi-state-restraints.refmac.params"

    for simul_occ in np.arange(start_occ, end_occ + step / 5, step):
        for starting_rand_occ in get_starting_occ(simul_occ, out_path, dataset_prefix):

            simulate_mtz = os.path.join(out_path,"{}_refine_occ_{}".format(dataset_prefix,str(simul_occ).replace(".", "_")),
                                        "{}_simul_{}.mtz".format(dataset_prefix, str(simul_occ).replace(".", "_")))

            out_path = os.path.join(out_path,"{}_refine_occ_{}".format(dataset_prefix,str(simul_occ).replace(".", "_")),
                                    "{}_expected_occ_{}_b_{}_supplied_occ_{}".format(dataset_prefix,
                                                                                     str(simul_occ).replace(".", "_"),
                                                                                     str(set_b).replace(".","_"),
                                                                                     str(starting_rand_occ).replace(".", "_")))
            if not os.path.exists(out_path):
                os.mkdir(out_path)
            os.chdir(out_path)

            refinement_random_occ_pdb = os.path.join(out_path,"{}_random_refine_occ_{}.pdb".format(
                                                                    dataset_prefix,
                                                                    str(starting_rand_occ).replace(".", "_")))

            sh_file = "{}_expected_occ_{}_b_{}_supplied_occ_{}.sh".format(dataset_prefix,
                                                                          str(simul_occ).replace(".", "_"),
                                                                          str(set_b).replace(".", "_"),
                                                                          str(starting_rand_occ).replace(".", "_"))

            out_prefix = "expected_occ_{}_supplied_occ_{}".format(str(simul_occ).replace(".", "_"),
                                                                  str(starting_rand_occ).replace(".", "_"))
            dir_prefix = "refine_" + out_prefix + "_cyc_20"

            print(out_path)

            os.system("cp multi-state-restraints.refmac.params multi-state-restraints-tmp.refmac.params")
            with open("multi-state-restraints-tmp.refmac.params",'a') as file:
                file.write('ncyc 20')

            with open(os.path.join(out_path, sh_file), 'w') as file:
                file.write("#!/bin/bash\n")
                file.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
                file.write(
                    "source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
                file.write("giant.quick_refine input.pdb={} input.mtz={} input.cif={} "
                           "input.params={} output.out_prefix={} output.dir_prefix={} ".format(
                    refinement_random_occ_pdb, simulate_mtz, input_cif, "multi-state-restraints-tmp.refmac.params",
                    out_prefix, dir_prefix,))

            os.system("qsub -o {} -e {} {}".format(
                os.path.join(out_path, "output_{}.txt".format(str(simul_occ).replace(".", "_"))),
                os.path.join(out_path, "error_{}.txt".format(str(simul_occ).replace(".", "_"))),
                os.path.join(out_path, sh_file)))

            out_path = os.path.dirname(os.path.dirname(out_path))