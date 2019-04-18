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


def submit_exhasutive_with_refmac_0(dataset_prefix, out_path, set_b=None):
    """ Qsub submission of sh files from occ_loop_merge_confs_simulate_with_refmac_0"""

    for lig_occupancy in np.arange(0.05, 0.96, 0.05):
        sh_file = "{}_occ_{}_b_{}.sh".format(dataset_prefix,
                                             str(lig_occupancy).replace(".", "_"),
                                             str(set_b).replace(".", "_"))

        os.system("qsub -o {} -e {} {}".format(
            os.path.join(out_path, "output_{}.txt".format(str(lig_occupancy).replace(".", "_"))),
            os.path.join(out_path, "error_{}.txt".format(str(lig_occupancy).replace(".", "_"))),
            os.path.join(out_path, sh_file)))


def refine_after_exhasutive_search(input_pdb, input_mtz, input_cif, refine_params, dataset_prefix, working_dir=None):
    os.chdir(working_dir)
    sh_file = "{}_quick_refine_exhaustive_search_minima.sh".format(dataset_prefix)
    out_prefix = "{}_refine_after_exhaustive_search".format(dataset_prefix)
    dir_prefix = out_prefix
    quick_refine_qsub(input_pdb=input_pdb, input_mtz=input_mtz, input_cif=input_cif, refine_params=refine_params,
                      sh_file=sh_file, out_prefix=out_prefix, dir_prefix=dir_prefix)


def quick_refine_qsub(input_pdb, input_mtz, input_cif, refine_params,
                      sh_file, out_prefix=None, dir_prefix=None):
    with open(os.path.join(out_path, sh_file), 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
        f.write(
            "source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
        f.write("giant.quick_refine input.pdb={} input.mtz={} input.cif={} "
                "input.params={} output.out_prefix={} output.dir_prefix={} ".format(input_pdb, input_mtz,
                                                                                    input_cif, refine_params,
                                                                                    out_prefix, dir_prefix))

    os.system("qsub -o {} -e {} {}".format(
        os.path.join(out_path, "output_{}.txt".format(str(out_prefix))),
        os.path.join(out_path, "error_{}.txt".format(str(out_prefix))),
        os.path.join(out_path, sh_file)))


def quick_refine_repeats(start_occ, end_occ, step, dataset_prefix, set_b, out_path, input_cif):
    params = "multi-state-restraints.refmac.params"

    for simul_occ in np.arange(start_occ, end_occ + step / 5, step):
        for starting_rand_occ in get_starting_occ(simul_occ, out_path, dataset_prefix):

            simulate_mtz = os.path.join(out_path,
                                        "{}_refine_occ_{}".format(dataset_prefix, str(simul_occ).replace(".", "_")),
                                        "{}_simul_{}.mtz".format(dataset_prefix, str(simul_occ).replace(".", "_")))

            out_path = os.path.join(out_path,
                                    "{}_refine_occ_{}".format(dataset_prefix, str(simul_occ).replace(".", "_")),
                                    "{}_expected_occ_{}_b_{}_supplied_occ_{}".format(dataset_prefix,
                                                                                     str(simul_occ).replace(".", "_"),
                                                                                     str(set_b).replace(".", "_"),
                                                                                     str(starting_rand_occ).replace(".",
                                                                                                                    "_")))
            if not os.path.exists(out_path):
                os.mkdir(out_path)
            os.chdir(out_path)

            refinement_random_occ_pdb = os.path.join(out_path, "{}_random_refine_occ_{}.pdb".format(
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
            with open("multi-state-restraints-tmp.refmac.params", 'a') as f:
                f.write('ncyc 20')

            with open(os.path.join(out_path, sh_file), 'w') as f:
                f.write("#!/bin/bash\n")
                f.write(
                    "export XChemExplorer_DIR=\"/dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer\"\n")
                f.write(
                    "source /dls/science/groups/i04-1/software/XChemExplorer_new/XChemExplorer/setup-scripts/pandda.setup-sh\n")
                f.write("giant.quick_refine input.pdb={} input.mtz={} input.cif={} "
                        "input.params={} output.out_prefix={} output.dir_prefix={} ".format(
                    refinement_random_occ_pdb, simulate_mtz, input_cif, "multi-state-restraints-tmp.refmac.params",
                    out_prefix, dir_prefix, ))

            os.system("qsub -o {} -e {} {}".format(
                os.path.join(out_path, "output_{}.txt".format(str(simul_occ).replace(".", "_"))),
                os.path.join(out_path, "error_{}.txt".format(str(simul_occ).replace(".", "_"))),
                os.path.join(out_path, sh_file)))

            out_path = os.path.dirname(os.path.dirname(out_path))


# Extra refinement for covalent ratios

start_xtal_num = 1905
end_xtal_num = 2005
in_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios"
out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios_phenix"
prefix = "NUDT7A-x"
qsub = False

if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    os.system('cp -a {}/. {}'.format(in_dir, out_dir))
#

xtals = []
for num in range(start_xtal_num, end_xtal_num + 1):
    xtal_name = prefix + "{0:0>4}".format(num)
    xtals.append(xtal_name)

for xtal_name in xtals:

    input_pdb = os.path.join(os.path.join(out_dir, xtal_name, "refine.pdb"))
    input_mtz = os.path.join(os.path.join(out_dir, xtal_name, "refine.mtz"))

    # if os.path.exists(os.path.join(out_dir, xtal_name,"refine_0002")):
    #     continue

    if not os.path.exists(os.path.join(out_dir, xtal_name, "dimple.pdb")):
        continue

    if not os.path.exists(os.path.join(out_dir, xtal_name, "multi-state-restraints.refmac.params")):
        continue

    print(xtal_name)

    # f = open(os.path.join(out_dir, xtal_name,
    #                  "multi-state-restraints.phenix.params"),"r")
    # lines=f.readlines()
    # f.close()
    # # while lines[-1].startswith("NCYC"):
    # #     lines.pop()
    # # if lines[-1].startswith("occupancy refine"):
    # #     lines.pop()
    # f = open(os.path.join(out_dir, xtal_name,
    #                  "multi-state-restraints.phenix.params"),"w")
    # f.writelines(lines)
    # f.write("strategy=occupancies")
    # f.close()

    refine_folders = [name for name in os.listdir(os.path.join(out_dir, xtal_name))
                      if os.path.isdir(name) and name.startswith('refine')]

    refine_fol_nums = [num[-4:] for num in refine_folders]
    final_refine = max([int(i) for i in refine_fol_nums])

    out_prefix = "refine_{0:0>4}".format(final_refine)

    cmds = "source /dls/science/groups/i04-1/software/" \
           "pandda-update/ccp4/ccp4-7.0/setup-scripts/ccp4.setup-sh \n"

    cmds += "cd {}\n".format(os.path.join(out_dir, xtal_name))

    # output.out_prefix =\"{}\"
    # cmds += "ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive/utils/quick_refine.py " \
    #         "{} {} {} params={} program=phenix\n".format(
    #     input_pdb,
    #     input_mtz,
    #     os.path.join(out_dir, xtal_name, "*.cif"),
    #     os.path.join(out_dir, xtal_name,
    #                  "multi-state-restraints.phenix.params"))
    #

    cmds += "ccp4-python /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive/utils/quick_refine.py " \
            "{} {} {} params={} program=phenix args=\"main.number_of_macro_cycles=1\" \n".format(
        input_pdb,
        input_mtz,
        os.path.join(out_dir, xtal_name, "*.cif"),
        os.path.join(out_dir, xtal_name,
                     "multi-state-restraints.phenix.params"))

    if qsub:
        f = open(
            os.path.join(out_dir,
                         xtal_name,
                         "{}_quick_refine.sh".format(xtal_name)), "w")

        f.write(cmds)
        f.close()

        os.system('qsub {}'.format(os.path.join(out_dir, xtal_name, "{}_quick_refine.sh".format(xtal_name))))
    else:
        print(cmds)
        os.system(cmds)

    exit()
