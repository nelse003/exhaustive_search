import os

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