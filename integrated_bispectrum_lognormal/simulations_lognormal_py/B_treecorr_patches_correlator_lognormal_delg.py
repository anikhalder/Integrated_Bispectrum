import time
start = time.time()

import numpy as np
import healpy as hp
import treecorr
import os
import sys

def createFolder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

if (sys.argv[1] == str(250)):
    ### Parameters to change according to patch size and count
    # Make 20 patches (discs) of 250 sq. degree pixels
    sq_degrees = 250
    #patch_radius = 0.155 #rad
    patch_count = 20
    #filepath = '../simulations_output/250_sq_degrees_20_patches/'
    maps_count = 10

elif (sys.argv[1] == str(50)):
    ### Parameters to change according to patch size and count
    # Make 100 patches (discs) of 50 sq. degree pixels
    sq_degrees = 50
    #patch_radius = 0.069 #rad
    patch_count = 100
    #filepath = '../simulations_output/50_sq_degrees_100_patches/'
    maps_count = 10

elif (sys.argv[1] == str(10)):
    ### Parameters to change according to patch size and count
    # Make 500 patches (discs) of 10 sq. degree pixels
    sq_degrees = 10
    #patch_radius = 0.031 #rad
    patch_count = 500
    #filepath = '../simulations_output/10_sq_degrees_500_patches/'
    maps_count = 10

elif (sys.argv[1] == str(5)):
    ### Parameters to change according to patch size and count
    # Make 1000 patches (discs) of 5 sq. degree pixels
    sq_degrees = 5
    #patch_radius = 0.022 #rad
    patch_count = 1000
    #filepath = '../simulations_output/5_sq_degrees_1000_patches/'
    maps_count = 10
    
else:
	raise Exception('Choose correct patch size!')

def calculate_patch_radius(patch_area_sq_degrees):
    return np.arccos(1-patch_area_sq_degrees*np.pi/(2*180*180))

def calculate_patch_area(patch_radius): 
    # patch is basically a spherical cap centered around the North Pole on a sphere of unit radius
    r = 1
    return 2*np.pi*r*r*(1-np.cos(patch_radius)) # from Wikipedia

patch_radius = calculate_patch_radius(float(sys.argv[1]))
filepath = '../simulations_output/'+str(sq_degrees)+'_sq_degrees_'+str(patch_count)+'_patches/'

################################################################

## TreeCorr
# for auto-correlation of the lognormal density fluctuations
# also multiply by mean density in patch to get integrated 3-pt function i_Xi

for j in range(maps_count):
    #print('Lognormal Map # ',str(j+1))
    filepath_map = filepath+'lognormal_map_'+str(j+1)+'/'
    createFolder(filepath_map+'B_treecorr_patches_correlated/')

    for i in range(patch_count):
        #print('Patch # '+str(i+1))
        density_fluctuations = hp.read_cl(filepath_map+'A_healpy_patches_produced/del_lognormal_patch_'+str(i+1)+'.fits')

        density_fluctuations_RA = density_fluctuations[0,:]
        density_fluctuations_dec = density_fluctuations[1,:]
        density_fluctuations_del_g = density_fluctuations[2,:]
        
        mean_del_g = np.mean(density_fluctuations_del_g)
        var_del_g = np.var(density_fluctuations)

        cat = treecorr.Catalog(ra=density_fluctuations_RA, dec=density_fluctuations_dec,
                                ra_units='rad', dec_units='rad', k=density_fluctuations_del_g)

        kk = treecorr.KKCorrelation(min_sep=1, max_sep=150, nbins=30, sep_units='arcmin')
        kk.process(cat) 
        theta_tc = kk.meanr # tc stands for treecorr
        xi_tc = kk.xi
        mean_del_g_vec = np.ones(theta_tc.size)*mean_del_g
        var_del_g_vec = np.ones(theta_tc.size)*var_del_g
        i_Xi = xi_tc*mean_del_g # integrated 3-pt correlation function
        
        dat = np.array([theta_tc, xi_tc, mean_del_g_vec, var_del_g_vec, i_Xi]) # the mean_del_g_vec, var_del_g_vec is the same value for all angles (just for output purposes putting into column)
        dat = dat.T
        np.savetxt(filepath_map+'B_treecorr_patches_correlated/correlation_lognormal_patch_'+str(i+1)+'.txt', dat, delimiter = ' ')


end = time.time()
print('\nTREECORR --> Time taken for execution (seconds): ', end - start)