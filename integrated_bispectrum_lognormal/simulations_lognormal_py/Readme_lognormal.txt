Instructions:

1. Produce lognormal maps using lognormal_maps_producer.sh - set the number of lognormal maps (maps_count) to be produced by flask using flask_lognormal.config file

2. Depending on the size patch you need in each map 100 patches of 50 sq degrees (each patch) or 20 patches of 250 sq degrees (each patch) - execute the following 3 python files:

A_healpy_patches_producer_lognormal_delg.py <sq_degrees_patch>
B_treecorr_patches_correlator_lognormal_delg.py <sq_degrees_patch>
C_patches_analyser_lognormal_delg.py <sq_degrees_patch>




### Parameters to change according to patch size and count
# Make 20 patches (discs) of 250 sq. degree pixels
sq_degrees = 250
patch_radius = 0.155 #rad
patch_count = 20
filepath = '../output/250_sq_degrees_20_patches/'
maps_count = 10


### Parameters to change according to patch size and count
# Make 100 patches (discs) of 50 sq. degree pixels
sq_degrees = 50
patch_radius = 0.069 #rad
patch_count = 100
filepath = '../output/50_sq_degrees_100_patches/'
maps_count = 10


### Parameters to change according to patch size and count
# Make 500 patches (discs) of 10 sq. degree pixels
sq_degrees = 10
patch_radius = 0.031 #rad
patch_count = 500
filepath = '../output/10_sq_degrees_500_patches/'
maps_count = 10
