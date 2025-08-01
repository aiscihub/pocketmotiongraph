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

load "data/ATRF_ASPFU_relaxed.pdb", protein
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

load "data/ATRF_ASPFU_relaxed.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [19757,11108,9113,9117,11005,9109,9184,9069,10958,10941,10896,21698,9256,9246,9247,9248,9250,21565,21620,21626,9131,9135,9124,19756,21607,21602,19740,19681,19685,21702,21605,19798,19802,21555,9188,9190] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.408,0.702]
select surf_pocket2, protein and id [21452,21502,21503,21507,21509,21511,21527,21531,21203,21523,21574,21578,21595,21266,21582,21975,21976,23893,23896,23897,23940,21305,21309,21313,8682,8793,8729,23944,8657,8661,8665,8733,23675,23679,21251,21253,23677,23596,22030,22029,21140,21142,22060,22063,21196,21199,23541,23520,23524,23464,23468] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.475,0.902]
select surf_pocket3, protein and id [9011,9015,9019,19559,19560,19563,11083,11087,11091,9003,9006,9009,9010,9013,9017,12055,19511,11142,11139,12071,19498,23135,23139,23057,23061,21738,23174,23172,23179,23196,19556,19508,21678,19612,19613,12355,12359,12280,12276,12335,12339,12357,12274,9058,9000,11591,11595,12338,12332] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.329,0.702]
select surf_pocket4, protein and id [9320,9339,9340,9341,9325,8616,8622,9293,8496,8502,8544,8548,8550,8443,8446,8452,8456,9844,9972,9904,10650,9863,9917,10570,9382,9386,9389,9385,8367,8368,8371,9789,4292,8437,8377,8427,8433,9795,4306,4317,4321,8430,8426,8441,9785,9792,9845,9841,9849,21378,21443,9318,21382,5141] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.361,0.373,0.902]
select surf_pocket5, protein and id [3506,6083,6085,6089,3489,3519,3521,3525,4219,4218,16014,3509,3554,3559,3565,3577,5558,3563,5551,5552,3543,3548] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.310,0.278,0.702]
select surf_pocket6, protein and id [19511,11142,11139,11135,19559,19563,11190,11186,11184,11246,11247,11243,19508,19515,19573,19574,19575,19485,19451,19526,19529,19368,19395,19403,11172] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.451,0.361,0.902]
select surf_pocket7, protein and id [12734,10526,12615,12679,9812,9879,9881,9883,10478,12563,12619,12606,12612,10459,10463,10477,10517,10521,10522,10524,12661,12658,12665,10472] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.392,0.278,0.702]
select surf_pocket8, protein and id [14453,14455,14457,14461,16901,5622,14497,14501,14513,16903,14472,14474,14479,14484,15101,16412,14527,5197,14506,15104] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.553,0.361,0.902]
select surf_pocket9, protein and id [16944,16962,16969,16544,16548,16551,16983,16989,16553,16557,16499,16496,16515,16518,16549,6985,6742,6745,16970,6904,6908,6989,7027,6982,6838,6842,6839,6843,6845] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.471,0.278,0.702]
select surf_pocket10, protein and id [3509,3513,3493,3512,3516,16032,6951,6949,13582,13586,13589,13590,6950,13554,16490,16521,16524,16030,16031,16519,16507,16506,16510,16042,15998,16046,3486,16117,15809,16564,16568,13597,13601,13576,13579,13592,6976] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.655,0.361,0.902]
select surf_pocket11, protein and id [10686,10664,10921,10987,10598,10989,10993,11368,11425,12698,11418,12647,11420,12716,12722,10918,10925,10970] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.549,0.278,0.702]
select surf_pocket12, protein and id [13043,13045,15848,15851,13219,13520,13524,2548,13276,13159,13209,13212,13152,13036,13039,13535,13537,13540,13543,13272,13000,2543,2545,12992,12997] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.757,0.361,0.902]
select surf_pocket13, protein and id [23636,23640,20482,20486,23547,20421,23551,23552,23558,21116,23482,23436,23498,21061,20441,20380] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.627,0.278,0.702]
select surf_pocket14, protein and id [2081,1938,1997,2000,2007,2075,1950,1988,1992,1993,4751,1934,1128,1049,1050,1994,4875,1884,1924,1927,1942] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.859,0.361,0.902]
select surf_pocket15, protein and id [5642,16901,5618,5622,14453,15100,15101,16412,16451,16903,5203,5205,15117,15104,5197,5208,5209,4243,5201,5210,5212] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.690]
select surf_pocket16, protein and id [12939,12938,9509,9510,9514,9515,9517,9583,9587,9591,9524,4170,4195,4199,3633,3636,3745,4166,4161,4164,3706,3708,3707] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.839]
select surf_pocket17, protein and id [4559,4587,8150,8114,4575,5448,2270,4124,2267,7947,7951,7907,2090,2110,5442,4659,4596,2116,2125,2128,2131,7893,7944,7950] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.612]
select surf_pocket18, protein and id [16928,17652,17656,17657,17699,17700,6804,6844,6845,5659,6098,16909,6092,6094,6095,16908,16924,16935,16946,16951,16955,16937,16912,16940,5663,5704,6144,17701,5714,6119,6117,5656] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.733]
select surf_pocket19, protein and id [1965,4891,5868,1936,4818,4834,2009,2097,2030] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.533]
select surf_pocket20, protein and id [2792,2685,2688,2692,2744,2878,10858,10877,2740,19114,19118,10764,10790,10797,10807,10806,10746,10760,10808,19051,19053,19056,19059,19045,2798,2612] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.361,0.631]
select surf_pocket21, protein and id [21392,21396,23960,23964,8533,23957,23958,8602,21333,21336,23889,23897,21309,21313] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.278,0.451]
select surf_pocket22, protein and id [11850,23010,23016,23031,23034,23035,23055,23061,11799,11881,22706,11804,12039,11811,12043,22816,22819] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.361,0.529]
select surf_pocket23, protein and id [3166,2498,3155,3197,2491,3521,3525,12975,12970,2534,13555,3577,3193,3760,3764,2486,3749,3187,3190] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.278,0.373]
select surf_pocket24, protein and id [18905,18909,20494,20504,20562,19099,2653,19095,20516,20506,20509,19103,20493,18991,18999,18900,18903] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.361,0.427]
select surf_pocket25, protein and id [9792,9497,9501,9693,9732,9749,9803,9742,9745,10635,12785,10703,12841,12838,9401,9448,9454,10707,9397,10698,10700,10636,10639,9748,9696] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.278,0.294]
select surf_pocket26, protein and id [21456,21266,21452,21531,21305,21309] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.392,0.361]
select surf_pocket27, protein and id [23675,23679,23867,23871,23677,23600,23655,23657,23659,23844] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.345,0.278]
select surf_pocket28, protein and id [19492,19498,23135,23057,23061,19502,23031,23035,12043,19563,22436,22437,12071,12039] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.902,0.498,0.361]
select surf_pocket29, protein and id [18631,20186,23783,20282,20297,20303,20316,18691,18695,18630,18644,18648,18741,18744] 
set surface_color,  pcol29, surf_pocket29 
set_color pcol30 = [0.702,0.427,0.278]
select surf_pocket30, protein and id [20168,20255,15359,20251,15347,16182,15362,18867,18875,15298,15295,15155,15170,15215,15219,15165,20154] 
set surface_color,  pcol30, surf_pocket30 
set_color pcol31 = [0.902,0.600,0.361]
select surf_pocket31, protein and id [19872,9315,9302,9303,9309,21428,21424,9321,19864,19868,21368,9370,9373,10777,19867,19927,10831,9320] 
set surface_color,  pcol31, surf_pocket31 
set_color pcol32 = [0.702,0.506,0.278]
select surf_pocket32, protein and id [9944,9946,10392,10384,10388,10383,10328,10308,10332,10324,10004,10390,10386] 
set surface_color,  pcol32, surf_pocket32 
set_color pcol33 = [0.902,0.702,0.361]
select surf_pocket33, protein and id [11163,12459,12461,11402,11406,12530,11500,12463] 
set surface_color,  pcol33, surf_pocket33 
set_color pcol34 = [0.702,0.584,0.278]
select surf_pocket34, protein and id [7639,7661,2278,2282,7776,3361,3363,7582,7586,3903,2312,7757,7763,7706] 
set surface_color,  pcol34, surf_pocket34 
set_color pcol35 = [0.902,0.804,0.361]
select surf_pocket35, protein and id [14234,14242,14213,14217,14239,14244,18254,18185,18182,18100,18106,17218,17223,17230,17226,18077,18083,14205] 
set surface_color,  pcol35, surf_pocket35 
set_color pcol36 = [0.702,0.663,0.278]
select surf_pocket36, protein and id [19869,19872,9309,21427,21490,9242,21487,9250,21549,21465] 
set surface_color,  pcol36, surf_pocket36 
set_color pcol37 = [0.894,0.902,0.361]
select surf_pocket37, protein and id [11316,11320,19316,19318,19251,19292,11349] 
set surface_color,  pcol37, surf_pocket37 
set_color pcol38 = [0.655,0.702,0.278]
select surf_pocket38, protein and id [17752,17746,14436,14439,5651,5652,5673,5653] 
set surface_color,  pcol38, surf_pocket38 
   

deselect

orient
