io_phil = libtbx.phil.parse("""
input{
    pdb = None
        .type = path
    mtz = None
        .type = path
    xtal_name = None
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
                            """, process_includes=True)