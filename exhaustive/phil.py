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
    log_name = "exhaustive_search"
        .type =str
}
settings{
    processes = 8
        .type = int
}
exhaustive{
    options{
        lower_occ = 0.0
            .type = float
        upper_occ = 1.01
            .type = float
        step = 0.05
            .type = float
        lower_u_iso = 0.2
            .type = float
        upper_u_iso = 1.21
            .type = float
        buffer = 0
            .type = float
        csv_name = 'u_iso_occupancy_vary'
            .type = str
        grid_spacing = 0.25
            .type = float
        generate_mtz = False
            .type = bool
    }
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
    options{
        set_b = None
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
    }
}
""", process_includes=True)

def prepare_validate_phil(master_phil):

    """ Add bound state and ground state paths if they do not exist"""

    params = master_phil.extract()

    if params.validate.input.bound_state_pdb_path is None:
        params.validate.input.bound_state_pdb_path = os.path.join(params.input.in_path,
                                                                  params.validate.input.bound_state_pdb_name)
    if params.validate.input.ground_state_pdb_path is None:
        params.validate.input.ground_state_pdb_path = os.path.join(params.input.in_path,
                                                                   params.validate.input.ground_state_pdb_name)

    modified_phil = master_phil.format(python_object = params)
    modified_phil.show()

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

