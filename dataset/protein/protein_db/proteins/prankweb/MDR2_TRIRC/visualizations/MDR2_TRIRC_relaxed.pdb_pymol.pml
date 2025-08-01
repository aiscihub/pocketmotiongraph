from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/MDR2_TRIRC_relaxed.pdb", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta




# SAS points

load "data/MDR2_TRIRC_relaxed.pdb_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [15526,15530,16084,16091,16092,16095,15481,1687,1691,5872,5815,5825,5831,5838,5812,5767,15478,2422,15546,15548,15550,13815,16104,13852,16049,16101,16105,12017,12021,5314,5315,5339,5876,5385,12075,13856,13795,11990,11985,12704,16046,12072,5774,16050,16051,16053,16055,16059,12133,12137,16057,15594,1737,1738,1740,1742,1746,5927,5869,5870,5874,5878,5865,5215,5219,5269,5268,16139,14387,15489,15435,15441,3500,15485,15486,3493,3551,3553,3556,3603,5830,5882,5885,5889,2563,2499,3542,3545,3549,2549,2495,13799,13798,12770,3490] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.369,0.702]
select surf_pocket2, protein and id [13742,16208,16135,5206,5210,5213,5217,5211,5215,5219,13746,13799,13688,12804,12807,12808,13693,12773,12812,12765,12866,5150,5093,5096,5143,5144,5147,5158,5160,12818,12814,5104,12859,12876,13686,13689] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.373,0.902]
select surf_pocket3, protein and id [13374,13452,13394,13375,13378,13379,13391,16449,16446,16455,13300,13302,13431,13436,13486,13490,13491,13493,13497,13499,13495,16459,13179,16374,16444,16371,16425,16424,16428,13554,16375,16376,16377] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.349,0.278,0.702]
select surf_pocket4, protein and id [14860,7194,7248,7346,7197,14842,14849,6711,6717,7366,7238,7855,9243,7830,14873,7859,7861,7866] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.553,0.361,0.902]
select surf_pocket5, protein and id [4582,4585,4593,17486,4575,4579,18113,17635,19392,18088,17531,17541,17639,18135,18124,17659,17539,17544,17489,4607] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.510,0.278,0.702]
select surf_pocket6, protein and id [13579,13583,13586,13589,13592,5084,12985,5101,12938,12939,12981,13641,12923,4330,4336,4312,4379,4383,4315,13067] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.757,0.361,0.902]
select surf_pocket7, protein and id [7189,7194,7198,7160,7162,7165,7167,7136,7133,7142,7144,7147,7148,6836,7212,6713,6717,7140,9257,7185] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.671,0.278,0.702]
select surf_pocket8, protein and id [11342,11348,16549,16557,13279,16561,17005,16960,17648,17653,17645,17664,17668,17672,16578,17669,16948,11413,11335,17001,17004,11367,13284,13253,13286,16518,13254,13248,13252,13283,13347,11291,11298,17705] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.839]
select surf_pocket9, protein and id [17425,17428,17434,17436,17439,17440,16992,16996,17452,17477,17504,17481,17454,17457,17459,19406,17486,17490] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.573]
select surf_pocket10, protein and id [2440,2442,2444,15457,2435,2436,2438,15517,15519,15515,15568,15570,2388] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.631]
select surf_pocket11, protein and id [11646,13675,11654,11656,11659,11701,11831,13657,13658,13668,13719,13723,11830] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.412]
select surf_pocket12, protein and id [14355,14359,14297,15417,14234,14291,14314,15513,15457,15517,15519,14240,15515,15410] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.427]
select surf_pocket13, protein and id [5217,5252,5269,5244,5869,5874,5878,5860,5854,4033,4034,5876,5872] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.306,0.278]
select surf_pocket14, protein and id [7966,15068,7945,7950,14968,14771,14727,14731,14735,14723,14762,14772,14777] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.498,0.361]
select surf_pocket15, protein and id [11856,11863,11926,11928,11930,11897,11902,11906,11910] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.467,0.278]
select surf_pocket16, protein and id [3372,14577,2684,3363,3376,3367,3378,3386,3439,14502,14505,2616,2619,15325,3423,14501,3430,15257] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.702,0.361]
select surf_pocket17, protein and id [15836,15839,12274,1893,15821,1965,1985,1964,1969,12209,15848,15853,12206,12208,12253,1907,12256,12257] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.624,0.278]
select surf_pocket18, protein and id [16472,16475,16479,11623,11627,11631,11687,11689,16468,13097,16493,11681,11685,13091,11680,11666,11669,11670,11671,16473] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.894,0.902,0.361]
select surf_pocket19, protein and id [20469,19823,20536,19836,19838,19840,17474,17408,20532,19833,19834] 
set surface_color,  pcol19, surf_pocket19 
   

deselect

orient
