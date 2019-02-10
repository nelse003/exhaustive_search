set normalize_ccp4_maps, 0

load /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios/NUDT7A-x1915/refine.pdb, refine-coot
select lig_cys, resi 73 or chain E
select surround_10, refine-coot within 10 of lig_cys
remove solvent
hide everything, refine-coot
color grey40, refine-coot
show lines, surround_10
show sticks, lig_cys
util.cbag lig_cys
bg_color white

load /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios/NUDT7A-x1915/coot_map.ccp4, coot_map
map_double coot_map, -1
isomesh map, coot_map, 0.139, surround_10

color skyblue, map

set mesh_width, 0.5
set ray_trace_fog, 0
set depth_cue, 0
set ray_shadows, off
set transparency, 0.5

# edge view
set_view (\
    -0.470849186,    0.065802455,   -0.879756033,\
    -0.594833195,    0.712762773,    0.371671110,\
     0.651516497,    0.698309839,   -0.296461374,\
     0.000064045,   -0.000038207,  -49.934093475,\
    39.262619019,  -43.407310486,   77.554275513,\
    12.287309647,   87.581398010,  -20.000000000 )

ray 2056,2056
png image.png


