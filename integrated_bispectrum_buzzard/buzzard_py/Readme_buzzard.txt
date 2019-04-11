Instructions:

1. The projected matter density map from Buzzard (only one map!) is already present in data folder 

2. Depending on the size patch you need in each map - execute the following 3 python files:

A_healpy_patches_producer_buzzard_del.py <sq_degrees_patch> 
B_treecorr_patches_correlator_buzzard_del.py <sq_degrees_patch>
C_patches_analyser_buzzard_del.py <sq_degrees_patch>




### Parameters to change according to patch size and count
# Make 20 patches (discs) of 250 sq. degree pixels
sq_degrees = 250
patch_radius = 0.155 #rad
patch_count = 20
filepath = '../buzzard_output/250_sq_degrees_20_patches/'


### Parameters to change according to patch size and count
# Make 100 patches (discs) of 50 sq. degree pixels
sq_degrees = 50
patch_radius = 0.069 #rad
patch_count = 100
filepath = '../buzzard_output/50_sq_degrees_100_patches/'


### Parameters to change according to patch size and count
# Make 500 patches (discs) of 10 sq. degree pixels
sq_degrees = 10
patch_radius = 0.031 #rad
patch_count = 500
filepath = '../buzzard_output/10_sq_degrees_500_patches/'


### Parameters to change according to patch size and count
# Make 1000 patches (discs) of 5 sq. degree pixels
sq_degrees = 5
patch_radius = 0.022 #rad
patch_count = 1000
filepath = '../buzzard_output/5_sq_degrees_1000_patches/'
