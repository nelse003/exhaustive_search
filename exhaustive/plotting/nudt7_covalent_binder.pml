set normalize_ccp4_maps, 0

load /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios/NUDT7A-x1958/refine.split.bound-state.pdb, refine-coot
select lig_cys, resi 73 or chain E
select surround_10, refine-coot within 10 of lig_cys
select surround_5, refine-coot within 5 of lig_cys
remove solvent
hide everything, refine-coot
color deepteal, refine-coot
show cartoon, refine-coot
color main set, lig_cys
show sticks, lig_cys
util.cbag lig_cys
bg_color white

load /dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/covalent_ratios/NUDT7A-x1958/coot_map.ccp4, coot_map
map_double coot_map, -1
isomesh map, coot_map, 0.139, lig_cys, carve=2.2

color skyblue, map

set mesh_width, 0.5
set ray_trace_fog, 0
set depth_cue, 0
set ray_shadows, off
set transparency, 0.5

set_view (\
    -0.759656191,   -0.215470746,   -0.613592684,\
    -0.479598641,    0.822840691,    0.304819256,\
     0.439211220,    0.525836349,   -0.728415370,\
     0.000064045,   -0.000038207,  -49.934093475,\
    39.262619019,  -43.407310486,   77.554275513,\
    12.287309647,   87.581398010,  -20.000000000 )

ray 2056,2056
png image.png