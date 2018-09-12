import libtbx.phil

copy_phil = libtbx.phil.parse("""
input{
    path = None
        .type = str
        .help = "Input path where input pdb & mtz files are stored" 
    base_pdb = None
        .type = str
        .help = "pdb (ground structure) from which atoms are copied from"
    pdb_style = "dimple.pdb"
        .type = str
        .help = "filename of repeated pdb" 
    mtz_style = "dimple.mtz"
        .type = str
        .help = "filename of repeated mtz"
    prefix = "NUDT7A-x"
        .type = str
        .help = "Prefix of xtal-name directory"
    start_xtal_number = None
        .type = int
        .help = "First xtal number "
    end_xtal_number = None
        .type = int
        .help = "Last xtal number"
    xtal_list = []
        .type = list
        .help = "Extra input xtal folders not covered by number range"
    atoms_new = []
        .type = list
        .help = "Residue id and chain of atoms to copy from base_pdb 
                to xtals in xtal list. in form ['E','1']"
    atoms_remove = []
        .type = list
        .help = "Residue id and chain of atoms to remvoe from xtals 
                in xtal list"
    link_record_list= []
        .type = list
        .help = "List of strings, formatted to PDB file 
                standard line format for a LINKR record.
                Use when copying a covalent ligand"
    cif = None
        .type = str
        .help = "Path to cif defining a compound added, for refinement"
}
output{
    out_dir = None
        .type = str
        .help = "Output directory"
    pdb = "dimple_with_lig.pdb"
        .type = str
        .help = "pdb filename with added and removed atoms"
    refine_pdb = "refine.pdb"
        .type = str
        .help = "pdb filename after refinement"
    multi_state_model_pdb = "multi-state-model.pdb"
        .type = str
        .help = "pdb filename after giant.merge_conformations"
}
settings{
    overwrite = False
        .type = bool
        .help = "Overwriting flag"
    qsub = True
        .type = bool
        .help = "Submit via qsub (at diamond)"
    ccp4_path = "/dls/science/groups/i04-1/software/" \
            "pandda-update/ccp4/ccp4-7.0/setup-scripts/ccp4.setup-sh \n"
        .type = str
        .help = "Path to source cpp4 setup"
    param_file = "multi-state-restraints.refmac.params"
        .type = str
        .help = "filename of parameter file for qucik refine"
}
""", process_includes=True)