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

load "data/MDR1_CRYNH_relaxed.pdb", protein
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

load "data/MDR1_CRYNH_relaxed.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [12952,12956,16893,16897,16901,660,6851,661,664,6855,13001,13003,13054,12986,12990,6252,17021,6327,6368,6374,6730,6798,13061,13065,6305,637,642,16958,646,6817,16962,2611,6810,6745,6749] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.376,0.702]
select surf_pocket2, protein and id [622,626,15228,15236,4497,4521,4562,6971,4566,6970,4622,16321,16325,629,3483,16367,16370,16371,16374,6909,6906,4643,4619,4633,4637,4673,2557,2560,2561,2564,6951,3627,2518,15322,15302,15306,17100,679,17047,17051,17055,723,726,17097,17104,698,6910] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.392,0.902]
select surf_pocket3, protein and id [734,738,739,14534,14538,14611,14593,737,761,14473,14472,14476,5346,5340,5342,5228,5290,5292,5233,6136,5296,13818,13821,13825,13827,14527,14531,14524,13815,14528,13761,13762,14596,5394,5326,6068,6071,6073,6137] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.329,0.278,0.702]
select surf_pocket4, protein and id [2192,2196,2141,2121,2125,2129,2145,2149,2430,4602,4604,2389,2372,2365,2373,2319,2323,2185,2189,2327] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.514,0.361,0.902]
select surf_pocket5, protein and id [742,710,5175,705,6910,6892,654,674,657,690,6234,6256,5225,5131,6233,5070,5074,5181,5188,5191,5192,5135,6938,6943,6944,6947,6954,6889,6896,6878,5091] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.471,0.278,0.702]
select surf_pocket6, protein and id [14167,14169,14172,14177,14183,17403,14164,14161,17398,17402,14126,12300,12301,12248,12302,12222,14165,14166,17335,17336,17339,12293,12297,14130,14132,17441,17406,17415,12254,14235,14236,14237,12315,14276,12345,12349] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.698,0.361,0.902]
select surf_pocket7, protein and id [15688,8354,15695,15701,8854,8858,8885,8829,8198,8242,8259,10260,15717,8865,8877] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.616,0.278,0.702]
select surf_pocket8, protein and id [18392,20356,18336,18338,18341,18342,18354,18356,18327,18379,18330,17910,17906,18361,18388,18383,18359,18406,18454] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.886,0.361,0.902]
select surf_pocket9, protein and id [13666,13658,13662,13678,13682,13684,13686,13609,13726,13728,12773,12779,13676,12783,13670,13776,6092,6099,12720,13656,13720,13725,13721,6163] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.639]
select surf_pocket10, protein and id [14442,14379,14384,14265,14261,14290,14340,14341,17286,14445,14059,14187,17360,13995,17285,17356,17351,17332,17336,17339] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.729]
select surf_pocket11, protein and id [15627,15629,8964,8965,8971,15635,15577,15583,15630,15580,15584,15586,8944,8949,15937,15941,9100,9102,15812,8974,8978,15870,15874,15570,15571] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.494]
select surf_pocket12, protein and id [8202,10274,8144,8146,8153,8189,7724,8148,8151,8152,8164,7847,8198,8169,8171,8166,8216,7851] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.545]
select surf_pocket13, protein and id [14665,14661,12502,14548,14550,14554,17173,17110,17126,17138,14624,14628,12498,14572,14576] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.353]
select surf_pocket14, protein and id [2116,2120,7085,2123,2127,7019,2125,2129,4533,4578,4602,6983] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.361]
select surf_pocket15, protein and id [4436,4437,4433,4501,3707,3759,15363,15425,4384,4386,15376,15377,4482,4486,4500,4502,4432,4441,4446,4443,4450] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.353,0.278]
select surf_pocket16, protein and id [734,14534,14538,14611,733,737,745,755,718,729,719,17116,716,17120,14596,14485,17161,14468,14479,17166,17207,17165,17167] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.545,0.361]
select surf_pocket17, protein and id [5158,5127,6111,6165,6167,6213,5154,5162,5164,6150,6154,6157,6217,6113,6117] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.494,0.278]
select surf_pocket18, protein and id [16391,16397,16346,16382,16339,16342,3496,3497,3501,3505,3492,16377,16444,3446,3452,16380] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.729,0.361]
select surf_pocket19, protein and id [14763,14722,17061,17064,17088,17090,15140,17068,14713,14706,14710,17030,14715,17034] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.639,0.278]
select surf_pocket20, protein and id [20340,20342,19006,19031,5625,5594,5601,18450,18454,18549,19042] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.886,0.902,0.361]
select surf_pocket21, protein and id [14698,17021,6333,12941,12895] 
set surface_color,  pcol21, surf_pocket21 
   

deselect

orient
