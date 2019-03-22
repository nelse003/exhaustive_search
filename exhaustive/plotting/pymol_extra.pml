extra commands

show surface, refine-coot


set normalize_ccp4_maps, 0

load /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios/NUDT7A-x1915/refine-coot.pdb, refine-coot
load /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios/NUDT7A-x1915/fo_fc_refine_coot.ccp4, coot_map_fofc
select lig_cys, resi 73 or chain E
hide everything, refine-coot
map_double coot_map_fofc, -1
isomesh map_fofc, coot_map_fofc, 0.127, lig_cys
color green, map_fofc
set mesh_negative_color, red
set mesh_negative_visible

bg_color white
set mesh_width, 0.5
set ray_trace_fog, 0
set depth_cue, 0
set ray_shadows, off
set transparency, 0.5
set_view (\
    -0.470849186,    0.065802455,   -0.879756033,\
    -0.594833195,    0.712762773,    0.371671110,\
     0.651516497,    0.698309839,   -0.296461374,\
     0.000064045,   -0.000038207,  -49.934093475,\
    39.262619019,  -43.407310486,   77.554275513,\
    12.287309647,   87.581398010,  -20.000000000 )

ray 2056,2056
png image_fofc.png