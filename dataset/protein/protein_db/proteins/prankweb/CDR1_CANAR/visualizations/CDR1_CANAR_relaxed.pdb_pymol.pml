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

load "data/CDR1_CANAR_relaxed.pdb", protein
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

load "data/CDR1_CANAR_relaxed.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [4133,4137,4171,7690,4175,4178,1748,1833,1835,4422,1681,7694,4206,1831,4207,4211,5010,5004,5007,5011,1611,1612,1695,4512,1617,4134,5015,5021,4538,5106,5549,5407,5463,5465,5468,5472,5524,5518,5454,5457,5461,4144,4152,5040,5045,5048,5049,5059,5073,1829,1986,5521,1965,1968,1978,5528,7499,7503,7506,7509,7512,1987,1929,4208,5109,5423,5403,5400,4475,4421,4465,4414,4516,4087,4153,4080,4083,5028,1351,1353,4155,4357,4359,4214] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.404,0.702]
select surf_pocket2, protein and id [12487,2780,2787,12482,3166,2783,3093,3094,3107,3109,3113,3148,3152,3125,3127,3132,3137,3143,3141,2213,12517,2758,3337,12480,3068,12501,3097,3101,3104,12492,12495,12510,5167,5170,5665,5669,5637,16195,12496,3224,3318,3348,3343,5625] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.459,0.902]
select surf_pocket3, protein and id [8853,8915,8919,8920,8923,8924,10213,8849,9414,9415,8850,8873,9416,8977,8978,8979,3894,3902,3928,8001,7983,7993,8090,8096,8164,8168,8100,8053,8059,8056,8015,8116,9411,9455,9458,9462,9530,9407,8002,8011,7982,7986,9475,7935,7930,7934,7938,9347,9351,7933,9410,8928,8935,8941,8985,9302,9059,9303,9305,9308,8903,21208,21210,21204,21255,8856,9055] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.310,0.702]
select surf_pocket4, protein and id [19518,21302,8790,21305,21377,8214,8215,8265,8755,8759,8248,21308,21311,21317,21306,8841,21250,21254,21255,21258,8854,8856,8149,8153,8157,8211,21253,21271,21321,8807,8814,9531,8789,8793,8805,8229,8216,8225,8164,8168,8147,8161,8167] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.380,0.361,0.902]
select surf_pocket5, protein and id [10710,19237,19280,19221,11573,11563,11575,22169,19208,19212,19209,19220,10607,19277,10615,8551,8555,8558,8546,8548,8554,11853,8534,8614,11773,8531,10658,11829,11580,11769,11833,11832,11826] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.341,0.278,0.702]
select surf_pocket6, protein and id [6573,15472,15470,15471,16276,6571,15720,15800,15808,16283,16287,16324,16331,15731,15735,15521,15687,3097,3101,3104,15719,15721,3081,3100,3084,13189,16242,3074,6536,6540,6541,16218,16229,6544,13134,13164,13185,6572,13128,13131,13188,13190,13234,3092] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.498,0.361,0.902]
select surf_pocket7, protein and id [21338,21407,21796,21713,21771,21799,21776,21714,21768,21027,23263,20971,23242,23192,23196,21327,21330,21334,21333,21386,21390,21394,21398,21072] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.431,0.278,0.702]
select surf_pocket8, protein and id [3810,3815,3823,3824,5167,5170,12500,15703,5665,5669,12496,12499,12457,3837,3224,3094,3107,3148,3152,12501,3097,3101,3805,3221,5156] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.616,0.361,0.902]
select surf_pocket9, protein and id [10703,10709,10738,10744,10706,10710,19230,19291,19226,22119,22123,19170,19172,19174,19240,19238,19179,19193,19198,19202,19132,19152,19163,19103,19108,19106,19114,10694] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.525,0.278,0.702]
select surf_pocket10, protein and id [531,535,1617,1622,1611,1612,4512,5334,4511,5454,5457,5461,1629,5403,5328,5398,5400,5379,4506,4503,558,575,579,1583] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.733,0.361,0.902]
select surf_pocket11, protein and id [9068,3290,3315,3320,9154,9157,3805,3832,3218,3221,3778,3782,3800,3318,9067,9078,9080,9151,3776,3287,3291,3295,12457,3837,3224,12441,12446,12449,12450,12453,3863,9021,9083,9087,9090] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.616,0.278,0.702]
select surf_pocket12, protein and id [20008,19822,19828,19901,19904,19903,19908,19838,23394,23490,23462,23480,23502,23505,23486,19980,20003,18367,18371,19957,18344,20013,20019,18347] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.851,0.361,0.902]
select surf_pocket13, protein and id [19325,22903,21533,21537,8551,8555,8556,8560,8587,8546,8514,8583,8566,19277,22849,22845,22893,21487] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.694]
select surf_pocket14, protein and id [16673,16679,16683,6316,6320,6380,6582,16206,16207,16210,16211,16259,16263,16264,16270,6579,6583,6586,6587,6588,16267,16274,6427,6431,6625,16656] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.831]
select surf_pocket15, protein and id [2902,2896,7086,7087,7095,7097,7140,2979,2980,2935,2938,6081,6078,2941,6044,6046,7030,7100,6058,6120,5845,6061,6063,6072,6050] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.600]
select surf_pocket16, protein and id [19517,19518,21302,8790,19521,21305,19510,8841,10473,19534,19463,8793] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.714]
select surf_pocket17, protein and id [16636,16640,16647,16649,16652,5700,6392,17457,17411,17458,5679,5683,5675,5687,5691,5695,5698,5262,6432,6433,6369,6351,6384] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.510]
select surf_pocket18, protein and id [19045,19049,10834,10836,10777,10910,10501,10450,10454,10849,10845,10853,10846,10519,10439] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.596]
select surf_pocket19, protein and id [14890,14966,14963,14833,15876,15880,14843,18610,14970,14823,14838,14809,19800,19893,19872,19875,19883,15014,15018,19897,19966,19970] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.416]
select surf_pocket20, protein and id [14498,14495,14497,14550,18196,18204,18209,14570,14580,14602,14581,14579,18277,18281,18212,14616,14670,14674,16013,14606,14609,14502,15985,13455,13452,15992,18215] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.361,0.478]
select surf_pocket21, protein and id [8690,21420,21373] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.278,0.325]
select surf_pocket22, protein and id [7677,7681,3705,3716,7655,12349,3630,3628,12394,12400,3634,3636,9077,9129,9132,9135,9126,9139,9141,9148,12437] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.361,0.361]
select surf_pocket23, protein and id [20378,20398,20402,20480,20332,20801,20802,20804,20748,20797,20806,20805,20808,20812,20814,20382,20741,20744,20745,20746] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.325,0.278]
select surf_pocket24, protein and id [14161,14168,14182,14149,14112,14116,14127,14134,14129,14108,14110,16615] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.478,0.361]
select surf_pocket25, protein and id [15016,18536,14954,15004,15006,15009,14955,14995,14998,18584,18590,18669,15057,15061,15213,15065,18527,18533,18540,18532] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.416,0.278]
select surf_pocket26, protein and id [9961,10029,10032,9979,9988,9965,9966,9967,9969,9973,9975,9488,9490,9499,9474,9540,9544,9470] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.596,0.361]
select surf_pocket27, protein and id [19045,19049,10834,19079,10777,19041,19094,19096,10853] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.510,0.278]
select surf_pocket28, protein and id [1639,1641,1651,515,519,1252,1243,1626,1717,1150,1198,1201,1711,1714,1629,1156,1162,1158,1655,1603,1644] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.902,0.714,0.361]
select surf_pocket29, protein and id [12513,12530,13119,2241,2270,2273,2276,12562,12566,12642,12646,13110,13114,13115,13117,13121,13125,2227] 
set surface_color,  pcol29, surf_pocket29 
set_color pcol30 = [0.702,0.600,0.278]
select surf_pocket30, protein and id [9364,12226,9380,9428,9432,10152,12151,10089,10093] 
set surface_color,  pcol30, surf_pocket30 
set_color pcol31 = [0.902,0.831,0.361]
select surf_pocket31, protein and id [6821,6797,13148,6818,6827,6815,13142,13146,13147,13190,6553,5950,5951,5936] 
set surface_color,  pcol31, surf_pocket31 
set_color pcol32 = [0.702,0.694,0.278]
select surf_pocket32, protein and id [9114,3879,4008,4058,3896,4007,4984] 
set surface_color,  pcol32, surf_pocket32 
set_color pcol33 = [0.851,0.902,0.361]
select surf_pocket33, protein and id [3467,3474,3250,2076,2071,3470] 
set surface_color,  pcol33, surf_pocket33 
   

deselect

orient
