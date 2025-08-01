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

load "data/CDR1_CANAR_auris_relaxed.pdb", protein
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

load "data/CDR1_CANAR_auris_relaxed.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5004,5010,5007,5011,4512,5407,1695,5524,5528,5468,5472,5463,5465,5518,5521,1986,5059,5109,5106,5549,1681,4471,4475,4516,4414,4418,4422,5028,4455,1831,1835,1935,1748,1353,1813,4220,1829,1833,1965,1978,1775,1968,1351,5400,5454,5457,4137,4178,7503,7506,7509,4171,4175,5073,1987,7512,4357,4359,1929,4208,4206,4211,4214,4216,5045,5015,5021,4134,5048,5049,4080,4083,4087,4534,4538,5423,5403] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.396,0.702]
select surf_pocket2, protein and id [7993,8001,8116,9415,8002,8011,8015,9411,9462,8985,9059,7924,7938,9302,9347,9351,7935,7934,7927,7930,7933,9400,8935,8941,8873,10213,9414,9416,9481,9458,9530,9475,9455,9407,9410,8096,8100,8056,8977,8978,8923,8979,8919,8920,8924,8928,8168,3928,7986,8053,3894,3902,7982,7983] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.443,0.902]
select surf_pocket3, protein and id [8149,8153,8157,8211,8214,8215,8854,21255,8856,21250,21254,21258,19518,21302,8790,21305,21377,8755,8759,8265,8841,8168,8229,8161,8167,8216,8228,8164,9531,8814,8225,8789,8805,8232,8807,8793,21253,21261,21306,21308,21271,21311,21317,21321] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.290,0.702]
select surf_pocket4, protein and id [3073,3077,3093,3094,3097,3101,3107,3109,5167,5665,5669,16199,16201,16195,3068,6481,5155,5170,5625,3137,3141,3143,3125,3130,3132,5637,12501,3152,3127,3148,3166,2783,3113] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.412,0.361,0.902]
select surf_pocket5, protein and id [8548,8551,8555,8558,8531,8534,8546,11573,10607,19237,11575,11580,10615,11853,8614,11773,11829,11769,11833,11832,19280,19277,19212,19223,19208,19209,19221,19220] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.373,0.278,0.702]
select surf_pocket6, protein and id [21327,21027,21334,21337,21338,21072,23263,23242,23246,21394,21389,21779,21713,21714,21768,21771,21776,21333,21407,21398,21400,21783,20967,20971,21796,21799,23192,23196,21386,21385] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.545,0.361,0.902]
select surf_pocket7, protein and id [16276,16280,16283,16287,16324,16331,15470,15472,15471,6571,6572,6573,16218,6540,6541,6544,16242,15719,15720,15721,15800,15808,15463,15466,13185,13188,13189,13131,13128,13134,13164,15687,15521,15731,15735,3074,16229,3097,3100,3104,3081,3084,3092] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.478,0.278,0.702]
select surf_pocket8, protein and id [3218,3224,3318,3348,3352,3343,2201,2780,2787,12495,12510,2758,12480,12482,12487,12492,12501,3152,3127,3148,3166,2783,12496,12499,3109,3143,3130,3337,12517,2276,2213] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.678,0.361,0.902]
select surf_pocket9, protein and id [12496,12499,12500,3824,3221,3810,3815,5156,12501,3152,3148,3094,3097,3101,3107,3823,5167,5665,5669,5170,3224,15703] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.580,0.278,0.702]
select surf_pocket10, protein and id [19163,19172,19240,19202,19238,19198,19226,22119,19114,19291,19170,19106,19132,19152,19174,19103,10703,10709,10694,10706,10710,19230,22123,10738,10744] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.816,0.361,0.902]
select surf_pocket11, protein and id [3863,9021,12446,12457,3218,3224,3318,12450,3856,3221,3805,3782,3800,3837,3832,3291,3295,9151,9154,9157,3776,3778,3287,3290,9067,9078,9080,9083,9087,9090,9181,12433,3315,9068,12441,12429,3320,12449] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.686,0.278,0.702]
select surf_pocket12, protein and id [23394,20008,23462,23480,23486,23490,20013,20019,19822,19828,19980,20003,18367,19838,18344,18347,19957,19903,19908,18371,19904,23502,23505] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.851]
select surf_pocket13, protein and id [14182,14116,14127,14134,14129,13814,13810,16123,14149,14139,14156,14161,14168,14143,14152,14110,14112,16615,14108] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.612]
select surf_pocket14, protein and id [19325,8556,8560,8566,19277,8587,8583,22903,21533,22889,22893,22849,22845,8551,8555,8546,8514] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.718]
select surf_pocket15, protein and id [16676,16679,16683,6625,6586,6588,16259,16263,16264,16270,6583,6579,6582,16267,16274,16210,6431,6427,16211,6380,16207,16656,6316,6320,16206] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.506]
select surf_pocket16, protein and id [21302,8790,8794,8797,21305,8841,10473,8793,19534,19463,21356,19510,19521] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.584]
select surf_pocket17, protein and id [12400,3630,3632,3634,3636,3716,9126,9129,9267,12349,9132,9135,3628,7655,7677,7681,9077,12394,12437,3705,3708,3295,9148,9139,9141] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.400]
select surf_pocket18, protein and id [5845,2979,6072,6078,7087,2896,2898,2902,6050,7030,7100,2980,7140,2935,2941,6046,6058,6061,6063,6044,7027,7032,6120] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.451]
select surf_pocket19, protein and id [14823,18610,15014,15018,19800,19893,19897,19966,14838,14843,14890,14963,14966,14970,15013,15009,19875,19883,19872,19970,14833,15876,15880] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.298]
select surf_pocket20, protein and id [5675,5679,5682,5683,5686,16636,16640,16647,6433,16649,16652,6432,16624,17457,6392,5262,5698,5700,5695,5691,5687] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.404,0.361]
select surf_pocket21, protein and id [20398,20480,20378,20382,20746,20748,20751,20741,20744,20745,20332,20402,20801,20802,20804,20805,20808,20797,20812,20814,20806] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.369,0.278]
select surf_pocket22, protein and id [9961,10029,10032,10089,9499,9966,9969,9965,9967,9973,9975,9488,9490,9540,9470,9474,9979,9988] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.537,0.361]
select surf_pocket23, protein and id [12209,10903,10910,10902,12205,10162,10519,10240,10454,10450,10501] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.471,0.278]
select surf_pocket24, protein and id [19045,19049,10775,10777,19096,19022,19086,19035,19041,19079,19094] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.675,0.361]
select surf_pocket25, protein and id [18277,18281,18209,18212,18215,18204,13452,14602,14670,14674,16013,13455,14498,14495,14502,15985,15992,14579,14580,14581,18196,14606,14609,14616] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.576,0.278]
select surf_pocket26, protein and id [15213,18583,18584,18590,18532,18527,18579,15004,15006,15009,15057,15061,14954,14995,14998,18533,18536,18540,15016] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.808,0.361]
select surf_pocket27, protein and id [12226,10152,12151,10089,10093,9324,9364,9380,9428,9432] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.678,0.278]
select surf_pocket28, protein and id [10511,10950,8780,10535,8733,10469] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.859,0.902,0.361]
select surf_pocket29, protein and id [3896,9114,3879,4008,4984,4007,4058] 
set surface_color,  pcol29, surf_pocket29 
   

deselect

orient
