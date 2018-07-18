import libtbx.phil
import os

master_phil = libtbx.phil.parse("""
input{
    pdb = None
        .type = path
    mtz = None
        .type = path
    xtal_name = None
        .type = str
    in_path = None
        .type = str
}
output{
    out_dir = "output"
        .type = str
    log_dir = "logs"
        .type = str
}
settings{
    processes = 8
        .type = int
    plot_dpi = 300
        .type = float
}
exhaustive{
    output{
        log_name = "exhaustive_search"
            .type =str
        csv_prefix = 'exhaustive_search'
            .type = str
        csv_name = None
            .type = str
    }
    options{
        lower_occ = 0.0
            .help = Lowest bound occupancy to check in the exhaustive search
            .type = float
        upper_occ = 1.00
            .help = Highest bound occupancy to check in the exhaustive search
            .type = float
        step = 0.05
            .type = float
            .help = Step size for u_iso and occupancy variation
        lower_u_iso = 0.2
            .type = float
        upper_u_iso = 1.20
            .type = float
        buffer = 0
            .type = float
        grid_spacing = 0.25
            .type = float
        generate_mtz = False
            .type = bool
        generate_map = False
            .type = bool
        convex_hull = True
            .type = bool
        convex_hull_buffer = True
            .type = bool
    }
}
select{
    resnames = DRG,FRG,LIG,UNK,UNL
        .help = 'Residues to generate constraint groups around for occupancy refinement (comma separated list of residue identifiers, i.e. resname=LIG or resname=LIG,UNL)'
        .type = str

    group_dist = 5
        .type = float
        .help = 'Distance to use when clustering atoms that should have the SAME occupancy'

    overlap_dist = 2
        .type = float
        .help = 'Distance to use when clustering atoms that should have occupancies that SUM TO LESS THAN ONE'

    exclude_altlocs = None
        .help = 'Exclude certain altlocs from occupancy groups (e.g. A or A,B)'
        .type = str

    complete_groups = True
        .help = 'Generate a set of fully constrained groups (that sum to unitary occupancy) when True. Generate a set of weaker constraints for overlapping atoms when False.'
        .type = bool

    verbose = True
        .type = bool
        
    coincident_cutoff = 0.05
        .help = 'RMSD Cutoff in Angstrom, for two structures considered coincident'
        .type = float
}
validate{
    input{
        bound_state_pdb_path = None
            .type = str
        bound_state_pdb_name = "refine.output.bound-state.pdb"
            .type = str
        ground_state_pdb_path = None
            .type = str
        ground_state_pdb_name = "refine.output.ground-state.pdb"
            .type = str
    }
    output{
        log_name = "validate"
            .type = str
        set_all_b_name_extension = "_set_all_b_"
            .type = str
        set_b_name_extension ="_set_b_"
            .type = str
    }
    options{
        set_b = None
            .type = float
        set_all_b = None
            .type = float
        step_simulation = 0.05
            .type = float
        start_simul_occ = 0.05
            .type = float
        end_simul_occ = 0.95
            .type = float
        buffer = 0
            .type = float
        overwrite = False
            .type = bool
        generate_ccp4 = False
            .type = bool
        use_qsub = True
            .type = bool
        qsub_error_prefix = "qsub_error_"
            .type = str
        qsub_out_prefix = "qsub_output_"
            .type = str
    }
}
repeat{
    input{
        compound_code = None
            .type = str
            .help = 'Compound code to refer to repeat soak compound. The one used in XCE/ sqlite database'
        database_path = None
            .type = path
            .help = 'Database path for sqlite databse from xce'
    }
}
""", process_includes=True)

def prepare_validate_phil(master_phil):

    """ Add bound state and ground state paths if they do not exist"""

    params = master_phil.extract()

    if params.validate.input.bound_state_pdb_path is None:
        print params.validate.input.bound_state_pdb_name
        print params.input.in_path
        params.validate.input.bound_state_pdb_path = os.path.join(params.input.in_path,
                                                                  params.validate.input.bound_state_pdb_name)
    if params.validate.input.ground_state_pdb_path is None:
        params.validate.input.ground_state_pdb_path = os.path.join(params.input.in_path,
                                                                   params.validate.input.ground_state_pdb_name)
    if params.validate.output.set_all_b_name_extension is None:
        params.validate.output.set_all_b_name_extension = str(params.validate.options.set_b).replace(".","_") +".pdb"
    else:
        params.validate.output.set_all_b_name_extension = params.validate.output.set_all_b_name_extension + \
                                                          str(params.validate.options.set_b).replace(".","_") +".pdb"

    if params.validate.output.set_b_name_extension is None:
        params.validate.output.set_b_name_extension = str(params.validate.options.set_b).replace(".", "_") + ".pdb"
    else:
        params.validate.output.set_b_name_extension = params.validate.output.set_b_name_extension + \
                                                     str(params.validate.options.set_b).replace(".", "_") + ".pdb"

    modified_phil = master_phil.format(python_object = params)

    return modified_phil

def check_input_files(params):

    if not os.path.exists(params.output.out_dir):
        os.mkdir(params.output.out_dir)

    assert params.input.pdb is not None, "Input pdb not supplied"
    assert os.path.exists(params.input.pdb), "Input pdb: \n{}\ndoes not exist".format(params.input.pdb)

    assert params.input.in_path is not None, "Input path not supplied"
    assert os.path.exists(params.input.in_path), "Input path:\n{}\n does not exist".format(params.input.in_path)

    assert params.input.mtz is not None, "Input mtz not supplied"
    assert os.path.exists(params.input.mtz), "Input mtz: \n{}\ndoes not exist".format(params.input.mtz)

