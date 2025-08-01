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

load "data/PDH1_CANGA_relaxed.pdb", protein
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

load "data/PDH1_CANGA_relaxed.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [1575,1630,1631,1634,1639,1586,1643,1645,1706,5021,1648,5478,4088,4286,4081,4085,4290,1702,1837,1798,1800,4116,4124,1686,4146,1804,4939,4957,4985,4953,4341,4345,4962,4944,1577,4951,4987,5426,5422,5481,5518,5521,1581,4049,7562,4991,7565,7568,1861,3620,1855,5018,4118,1819,1834,1847,4982] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.404,0.702]
select surf_pocket2, protein and id [12717,12710,12712,12716,2642,3033,2977,2991,12679,12690,2646,2038,2639,2597,2601,2605,2051,3198,5122,3007,3012,3015,3019,5119,5621,5623,5627,3005,2941,2971,2973,3001,12725,12726,12729,12731,12745,2958,12730,3183,3209,3197,3204,3085,3091,5583] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.459,0.902]
select surf_pocket3, protein and id [8667,8673,8674,8675,8677,8679,8685,10780,19787,19788,19784,19791,8665,8670,8627,8618,8621,11767,11771,11764,8702,8704,19745,19748,22664,11785,22657,11773,10772,10821,19836,23425,23370,23408,23415,23366,19718,19724,21976,21984,22028,22032,23410,8732,12067,12071,8650,8653,11980,11982,8657,12045,12049,12048] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.310,0.702]
select surf_pocket4, protein and id [21770,21772,8978,21806,8317,21823,8880,8874,8878,8253,8246,8238,8242,21763,8240,8232,8299,8302,8926,8914,8963,20095,8910,8911,8916,21810,21880,21876,21766,21812,21817,21814,21827,20026,20021,8256,8257,8260,8934,8928,8931,8314,8870,8252] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.380,0.361,0.902]
select surf_pocket5, protein and id [21586,21584,23795,23799,23811,22278,23791,21484,21488,22305,22308,23736,22277,22288,23740,22274,22286,21833,21836,21840,21844,21887,21895,22217,22221,22220,21538,21850,21839,21886,21890,21899,21901,21904,22216,22291] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.341,0.278,0.702]
select surf_pocket6, protein and id [16210,16212,16714,16724,6567,2956,6534,6538,6541,2961,2965,2938,2942,2945,2968,16711,2964,16178,13469,13472,13475,13457,13458,13459,5894,16012,16222,15951,16211,16299,16765,13591,6549,13595,15945,13518,13524,13528,13587,16226] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.498,0.361,0.902]
select surf_pocket7, protein and id [8148,8152,8199,8203,8090,8096,8103,9568,8073,8076,8087,9105,3815,9103,8079,9049,8256,8257,8260,8986,8995,9096,9045,8016,9464,9521,9524,9525,9510,9457,9460,9526,21691,21714,21721,21723,21736,8975,21722,8971,8180,8183,21774,21695,3819,3821,8187,24461,8173] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.431,0.278,0.702]
select surf_pocket8, protein and id [8862,8911,8918,21880,21876,8804,8808,8852,8858,8922,8874,8878,21924] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.616,0.361,0.902]
select surf_pocket9, protein and id [9439,9443,12459,12461,12372,9441,9474,12435,12437,12439,12441,12432,12433,10215,10219,10220,10222,10282,12313,12356,12312,12355,12360,12366,12428,9542,10224,10228,10226,9490] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.525,0.278,0.702]
select surf_pocket10, protein and id [19464,19462,19466,19530,19534,19520,10952,11026,10595,10599,10670,11007,11009,11090,11022,10617,10623] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.733,0.361,0.902]
select surf_pocket11, protein and id [3694,3664,3668,3689,3088,3091,3751,9147,9192,9196,9204,9207,9211,9277,9273,3662,3155,12658,12664,3180,3184,12667,9201,3185,9214,3723,9299,3719,12683,12684,12687,12648] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.616,0.278,0.702]
select surf_pocket12, protein and id [5108,3012,3015,3019,5119,5623,5627,2971,3711,3712,3713,2958,16194,12730,12726,12729,12731] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.851,0.361,0.902]
select surf_pocket13, protein and id [9413,9414,9066,9185,9178,9181,9464,9460,9048,9053,9059,9111,9056,9118,9121,9125,12559,12566,10387,12611,9372,12565,10391,10384,10390,9128,9131,9139,9375,9417,9420,9415] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.694]
select surf_pocket14, protein and id [10404,10371,10619,10293,10359,10617,10623,10690,11090,10353,12417,12471,10296,12421,11082] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.831]
select surf_pocket15, protein and id [16502,14955,14956,14969,13944,14974,18780,15056,15057,15058,15082,18700,18697,18692,18688,15096,16496,16498,14953,14960,18762,16470,16477,15144,15148] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.600]
select surf_pocket16, protein and id [4435,5329,5287,5336,4436,241,243,271,273,4430,5285,5286,5353,1510,1523,1525,1528,1491,1517,1521,247,5411] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.714]
select surf_pocket17, protein and id [14664,14650,14592,14594,14614,14616,14638,16605,17084,14580,14631,14634,14590,14611,14609] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.510]
select surf_pocket18, protein and id [13509,13425,2099,13044,12800,13508,13424,13435,13445,12757,2103,12754,12760] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.596]
select surf_pocket19, protein and id [5649,5653,5656,5658,5661,6371,17116,6420,6421,17121,5637,17928,17929,17882,17090,5214,17118,17887,17127,17093,17105] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.416]
select surf_pocket20, protein and id [22338,22239,22245,22255,22261,22358,22354,23614,22408,23441,22414,23542] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.361,0.478]
select surf_pocket21, protein and id [20479,19076,20388,20392,15505,15509,19065,15443,15446,15449,15450,16367,16371,15292,15307,20370,20305,15312,20378,15302,15367,20475,15504] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.278,0.325]
select surf_pocket22, protein and id [16692,17150,16741,16745,16746,17149,17142,6608,6616,16688,16689,6359,6415,6419,6420,6421,16693,17121,17125,6418,6556,6559,6563,6568,6562,6473,16749,16752,16756] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.361,0.361]
select surf_pocket23, protein and id [18816,18817,18834,18836,20466,20510,20517,20528,20522,20551,20554,20403,20406,24050,24031] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.325,0.278]
select surf_pocket24, protein and id [21460,23711,23757,23761,20731,23829,20727,21464,20795,20802,21393,21403,21406,21340,23649,21347,23707] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.478,0.361]
select surf_pocket25, protein and id [3388,3391,7267,7332,7328,7261,7271,7307,3387,3380,3383,2809,2812,2781,2788,2791,2794,2752,1900,2768,1903,1905] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.416,0.278]
select surf_pocket26, protein and id [10623,10684,10690,11090,10296,12421,11082,10670] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.596,0.361]
select surf_pocket27, protein and id [10944,10952,10935,10939,10940,10942,10943,10946,19580,19582,19573,19576,19569,19574,19578] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.510,0.278]
select surf_pocket28, protein and id [19578,19628,19632,19634,19636,19574,19582,10940,10946] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.902,0.714,0.361]
select surf_pocket29, protein and id [2195,19193,20769,20771,19185,19372,20814,20818,20761,20763,20738,20739,19106,19110] 
set surface_color,  pcol29, surf_pocket29 
set_color pcol30 = [0.702,0.600,0.278]
select surf_pocket30, protein and id [9443,12457,12461,7808,7828,9441,12511,12437] 
set surface_color,  pcol30, surf_pocket30 
set_color pcol31 = [0.902,0.831,0.361]
select surf_pocket31, protein and id [3599,3502,12581,7746,7750,7724,7727,9365,9256,9264] 
set surface_color,  pcol31, surf_pocket31 
set_color pcol32 = [0.702,0.694,0.278]
select surf_pocket32, protein and id [24459,8236,24450,24453,24457,8242,8240,21783,21790] 
set surface_color,  pcol32, surf_pocket32 
set_color pcol33 = [0.851,0.902,0.361]
select surf_pocket33, protein and id [12185,9987,9995,12191,11348,12130,11354,12133,9981] 
set surface_color,  pcol33, surf_pocket33 
   

deselect

orient
