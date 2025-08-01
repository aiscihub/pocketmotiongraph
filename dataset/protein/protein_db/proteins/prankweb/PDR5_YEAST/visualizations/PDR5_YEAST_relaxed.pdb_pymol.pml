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

load "data/PDR5_YEAST_relaxed.pdb", protein
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

load "data/PDR5_YEAST_relaxed.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [8169,8173,8747,8751,8745,8749,8205,21444,21448,21456,21462,21521,21517,21402,21405,21406,21410,21445,21449,8853,8847,8849,8111,8170,8799,8802,21566,19557,19553,19618,19622,8793,10551,8675,8679,19571,8741,8729,10547,8663,21499,21503,19670,19674,19675,19687,19690,19744,19617,8834,10436,10490,8789,8795,8781,8782,8785,19143,10499,19572,8086,8155,8158,8159,8160,8162,8110,8096,21466,8092] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.396,0.702]
select surf_pocket2, protein and id [4371,4375,4855,4851,4854,5258,5330,5331,5334,4026,4056,1614,1763,1767,5381,5385,1609,1620,1669,1673,5378,1800,5255,5262,5327,5318,5409,7407,4893,4917,4920,4923,7413,1820,1824,7410,4031,1797,1782,1761,4884,4889,4892,4309,4259,4255,3995,3998,4034,4198,4202,3956,3959,4864,3991,5364,1547,1539,1561,1606,1610,1650,4349,1499,1503,5309,5312,1556,1605,1618,4304,4306,4296,4346,4353] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.443,0.902]
select surf_pocket3, protein and id [2964,5010,2957,2960,3036,12594,3030,2886,2934,2916,2918,2941,2946,2950,12598,12597,12599,2882,2902,2903,2906,12564,12578,12580,12583,12585,12558,12590,12593,2922,13059,2603,12547,2978,3128,3154,3158,3143,12561,12565,12608,12615,2565,2600,2607,3139,3142,2569,2002,2007,2016,13040,3629,5021,5024,5515,2877,5505,5483,3620,3628,2871,2843,5509] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.290,0.702]
select surf_pocket4, protein and id [7912,8035,7943,8057,7990,9381,9392,21356,21363,21365,7995,8037,8039,8041,8043,3728,8976,7913,7916,7919,7927,3734,8911,8913,8916,8846,8857,8859,8920,8842,21360,8866,8870,9395,9396,9397,8107,8103,8053,9517,9440,9436,8805,8843,8111,8801,8802,8806,9513,9451,9457,9444,10055,8808,8810,8817] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.412,0.361,0.902]
select surf_pocket5, protein and id [19490,8547,8551,19442,19445,23084,23024,23028,23067,22927,11615,11620,11624,11630,11638,19402,11922,11928,10669,11896,11899,11902,11842,11840,11906,19388,19391,19372,19375,19386,22319,19378,22301,22315,11650,8549,21674,8544,8545,8572,8574,23066,23070,21626,8555,8535,8520,8483,8492,8498,10628,8602,8537,8540,8523] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.373,0.278,0.702]
select surf_pocket6, protein and id [9310,9314,9345,9368,9413,12299,12305,12307,12303,10138,10136,12224,12235,12245,12298,12301,12327,10072,12182,12225,10023,10033,12181,10073,10077,10081,10075,10079] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.545,0.361,0.902]
select surf_pocket7, protein and id [10880,10892,10899,10817,10825,10882,19179,10443,10447,10458,10518,10895,19195,19199,10963,10538,10955,19121,19185,19188,19127,10465,10471,19125,19129,19131] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.478,0.278,0.702]
select surf_pocket8, protein and id [8937,9283,9284,9056,9052,9285,8930,8982,8927,8989,8992,8996,9282,9286,9288,9291,9246,9331,9335,10235,12433,10245,8999,9002,9010,12479,9243,12426,12432,10244,10238] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.678,0.361,0.902]
select surf_pocket9, protein and id [9018,9062,9073,9075,9063,9148,9144,12551,12555,3578,3580,3584,3610,3097,3101,3105,3036,3030,3033,12526,12529,12532,3100,3125,3130,9078,9082,9085,9164,9170,3639,3664,3605] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.580,0.278,0.702]
select surf_pocket10, protein and id [13060,5785,5786,13071,13112,13070,13074,5787,13171,13128,13031,2567,2571,2573,2583,2645,2648,5811,5784,2599] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.816,0.361,0.902]
select surf_pocket11, protein and id [21485,21179,21481,21549,21866,21923,21922,21950,21953,21224,23427,23431,23450,23377,23381] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.686,0.278,0.702]
select surf_pocket12, protein and id [14257,14260,14284,14289,14296,14310,14229,14238,16774,14226,14240,14255,14277,14280,14262,14271,14244,14245,13944,16289,14236] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.851]
select surf_pocket13, protein and id [5528,5529,5532,16806,16808,16783,16795,6268,16811,6227,5525,6267,16799,17557,17562,5521,5533,5537,5541,5544,5546,5555,5116,17603,17604] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.612]
select surf_pocket14, protein and id [2906,2909,15892,15858,15902,15906,13111,15696,15627,15890,2910,16408,16392,15891,15979,16453,16449,15971] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.718]
select surf_pocket15, protein and id [12337,10467,10150,12263,12283,12287,12267,10258,10207,10225,12416] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.506]
select surf_pocket16, protein and id [1431,1468,1427,1417,1413,1402,436,413,4411,441,4378,4381,4340,4333,448,4336,1473,452,4413,4425,4407] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.584]
select surf_pocket17, protein and id [6265,6215,6262,6161,16376,16373,16815,16832,16842,6404,6327,6401,6407,6408,6411,6413,16433,16440,16430,16436,16838] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.400]
select surf_pocket18, protein and id [3185,3187,3181,3355,3359,3232,3231,3118,3429,3377,3378,3368,3426,3083,3086,3090,3361,3363,3366,3265,3268,3275,3267,3260] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.451]
select surf_pocket19, protein and id [8599,8404,8459,8463,11295,9742,9758,9760,9763,11255,11251,9737,9734,9746,8415,9782] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.298]
select surf_pocket20, protein and id [20376,20380,23463,23348,23352,23344,23399,23403,21103,23467,21107] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.404,0.361]
select surf_pocket21, protein and id [20158,20168,20174,23670,20049,20135,20052,20055,23689,23685,20112,18483,18484,18485,18502,18504] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.369,0.278]
select surf_pocket22, protein and id [9127,9135,12448,9236,3419,3515,3441,3435,12498,7595,7591,7569,7572,3504] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.537,0.361]
select surf_pocket23, protein and id [4349,456,5183,5189,1503,490,494,496,492,4348,476,477,554,1471,1485,1488,1492,5253,5255,5258] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.471,0.278]
select surf_pocket24, protein and id [6853,6938,6942,6857,2659,6784,6787,2694,1948,2656,6831,6837,5832] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.675,0.361]
select surf_pocket25, protein and id [18937,19824,18815,18893,18897,18930,19828,18918,18924,19756,19878,19882,19821,19812] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.576,0.278]
select surf_pocket26, protein and id [18987,2146,2167,18985,2147,19031,19034,19038,20470,20420,18862] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.808,0.361]
select surf_pocket27, protein and id [11403,11509,11516,11464,11360,11394,11406,11355,11359,11362,11366,11580,11364,11370,11727,11730,11536,11576,11533] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.678,0.278]
select surf_pocket28, protein and id [22065,21884,23102,21900,21903,21906,22009,22003,23255,23183] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.859,0.902,0.361]
select surf_pocket29, protein and id [1545,4237,1188,4190,4194,4198] 
set surface_color,  pcol29, surf_pocket29 
   

deselect

orient
