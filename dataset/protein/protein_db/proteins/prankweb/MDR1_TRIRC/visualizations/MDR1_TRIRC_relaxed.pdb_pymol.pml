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

load "data/MDR1_TRIRC_relaxed.pdb", protein
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

load "data/MDR1_TRIRC_relaxed.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5376,5384,5385,5341,5344,5348,5351,5379,5381,5438,5395,5471,5748,5907,5799,5846,5866,5873,4426,5357,4508,2190,2191,2257,2273,4463,4467,4492,2289,2296,5850,2299,2302,2310,4425,4429,4445,4470,4801,1991,2058,1987,2062,2065,5795,5844,5781,5364,4724,4730,4734,4675,2121,2125,2129,2143,4673,4542,4867,5727,5732,4044,2321,2311,2316,2320,2319,8111,8119,7865,7918,7866] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.400,0.702]
select surf_pocket2, protein and id [9201,9142,9191,9197,9088,9092,19620,10793,9244,9248,10891,10887,10842,10910,10948,10951,10955,19616,10786,10831,10836,9194,19677,19680,19684,21474,21528,21531,19745,19730,19734,19679,9214,9217,9218,9151,9161,8673,8677,9154,8617,8624,9260,21422,21429,21431,21441,21482,21470,8680,8681,21551,19795,21423,21427,19794,19798,8613,8608,8612,8615] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.451,0.902]
select surf_pocket3, protein and id [1315,1321,1323,1380,1327,1331,2103,1836,1842,1782,1784,4557,1771,1774,1778,1304,1305,1788,1787,4659,4669,4673,2111,2113,2116,2123,2127,4584,4580,4535,4538,2174,4714,4715,4760,4771,4717,2052,1839,2046,1884,1859,1860,1861,1855,4704,4707,4712] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.302,0.702]
select surf_pocket4, protein and id [8955,8958,11920,19451,22289,22293,22282,11910,22279,11907,22888,11028,19521,11082,19455,19457,19500,19496,11085,19460,11036,19568,8962,8965,8967,12107,12162,12103,11929,22276,12178,12184,11032,8953,9007,8938,8907,11916,11924,11512,19507,22979,22983,22987,22905,19443,19450,22884,22903,21601,21661,21664] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.396,0.361,0.902]
select surf_pocket5, protein and id [9330,9336,10604,9275,9279,9280,9793,9791,9792,9861,9380,9685,9689,8407,9384,9387,9268,9316,9320,9680,9684,9693,9730,9732,9735,9736,8463,8469,8389,8388,8392,8410,8423,8507,8508,8512,8514,8516,8522,8526,4233,8417,8408,8338,8341,8345,8419,9842,9838,8399] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.357,0.278,0.702]
select surf_pocket6, protein and id [3361,3374,3376,5497,5983,5987,3380,3394,3344,3392,4125,4126,4127,12818,12821,12822,12823,15962,15983,3364,5549,5486,3410,3415,3419,3397,3399,3404,3489,3492,4118,3056,3433,3052,12796,12801,12809] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.522,0.361,0.902]
select surf_pocket7, protein and id [11278,11281,11344,11266,11268,10868,10871,10872,10939,10851,10919,11285,19265,10859,10861,10864,10815,10553,10557,11336,10549,12531,12535] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.455,0.278,0.702]
select surf_pocket8, protein and id [4133,4136,4165,4171,3489,3492,4081,4085,4091,4103,4108,9537,9542,4079,9560,9564,9469,9419,9477,9479,9482,9486,9489,4139,12761,12765,3582,3548,3557,3551,12779,12782,12786] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.647,0.361,0.902]
select surf_pocket9, protein and id [21219,21456,23376,23380,21515,21920,21906,21917,21934,21937,21907,21175,21518,21120,21171,23316,23320,23395,21492,21506,21847,21498,21448,21452] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.549,0.278,0.702]
select surf_pocket10, protein and id [1771,1778,2155,2161,1732,2106,2107,1207,1213,1209,1272,1143,1073,1077,1078,1080,1082,1139,2174,2215,2220] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.773,0.361,0.902]
select surf_pocket11, protein and id [6734,6934,6937,6940,6941,6944,6946,6737,16471,16475,16470,16474,16477,6785,6781,16924,16930,16934,6677,16906,16910,16947,16953,6787,16521,16525,16528,16530,16529,16532,16536,16538] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.647,0.278,0.702]
select surf_pocket12, protein and id [9217,8565,9218,9257,9260,9264,21379,21420,21423,21427,21381,21384,21385,21429,21431,21441,21341,8551,8556,8557,8559,8569,8609,8613] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.898,0.361,0.902]
select surf_pocket13, protein and id [12669,12672,12675,9639,9643,10648,12676,10650,10653,10657,9397,9400,9403,9411,10642,10647,12720,12741,9450,9453,9376,9457,9685,9387,9336,9329,9390,9326,9680,9691,9693] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.655]
select surf_pocket14, protein and id [14358,14354,14398,14347,14356,14361,14407,14414,14375,14428,14395,14402,14373,14380,14362,16387] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.776]
select surf_pocket15, protein and id [14098,13821,18201,18129,18197,18050,14134,17224,14101,14126,14131,17220,18034,18180,18038,18051,18057,18044,17217,18028,14089] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.557]
select surf_pocket16, protein and id [12509,12511,12513,12515,12502,12506,12507,12462,12458,12531,12535,11334,11336,11340,11344] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.651]
select surf_pocket17, protein and id [11298,11404,11304,11318,11322,12211,12172,12193,12196,12343,11469,11470,11471,12272,12274,11415,11419,11425,11300,11090,11094] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.459]
select surf_pocket18, protein and id [18869,18870,18874,20473,20477,20483,20453,18868,18865,18949,18957,2682] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.525]
select surf_pocket19, protein and id [15084,15069,15072,15079,15220,16155,20086,15089,15137,15140,15144,15273,15278,15217,16151,18843,20008,18832,20004,20100] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.361]
select surf_pocket20, protein and id [10277,10336,10340,10281,9886,9947,10261,9884,10357,10351,10275] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.361,0.400]
select surf_pocket21, protein and id [4217,4220,4199,4233,8396,4239,8463,8469,8389,8388,8392,8393,9365,9366,9370,9380,9441,9384] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.298,0.278]
select surf_pocket22, protein and id [19548,19507,22979,22983,22987,22989,22993,23006,21601,21661,23029,23033,8967] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.443,0.361]
select surf_pocket23, protein and id [10851,11285,19261,19265,11218,11222,11266,19208,19254] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.392,0.278]
select surf_pocket24, protein and id [5528,16428,16430,16433,16434,16419,5972,5979,16443,5541,5535,5550,5565,6004,5552,16426,5519] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.569,0.361]
select surf_pocket25, protein and id [8254,8257,8263,8239,8242,8304,8245,8308,8380,4361,4363,4406,4611,4625,4628,4632,8252,8188,8235,4614,4354] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.490,0.278]
select surf_pocket26, protein and id [23411,23286,23358,20445,23338,21096,23294,20441] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.694,0.361]
select surf_pocket27, protein and id [14328,16883,17818,17749,17753,14350,14345,14331,5581,16867,16869,16863,5585,5586,17822,17815] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.588,0.278]
select surf_pocket28, protein and id [9633,7980,8003,8055,8059,8015,8162,8127,7997,8150,8141,8210] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.902,0.820,0.361]
select surf_pocket29, protein and id [22157,22161,22164,22099,20756,20754,20797,20793,20799,20815,22167,22175,19425] 
set surface_color,  pcol29, surf_pocket29 
set_color pcol30 = [0.702,0.686,0.278]
select surf_pocket30, protein and id [3158,6381,6385,7475,6171,6174,6379,6393,7517,7521,7522,891,6406] 
set surface_color,  pcol30, surf_pocket30 
set_color pcol31 = [0.855,0.902,0.361]
select surf_pocket31, protein and id [6014,6787,17698,17657,6705,6746,16906,16894,16901,16903,16912,16878,5993,6005,6009,6015,6017] 
set surface_color,  pcol31, surf_pocket31 
   

deselect

orient
