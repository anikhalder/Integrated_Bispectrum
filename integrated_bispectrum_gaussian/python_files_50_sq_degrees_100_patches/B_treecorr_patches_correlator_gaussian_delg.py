import time
start = time.time()

import numpy as np
import healpy as hp
import treecorr
import os

def createFolder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

### Parameters to change according to patch size and count
# Make 100 patches (discs) of 50 sq. degree pixels
sq_degrees = 50
patch_radius = 0.069 #rad
patch_count = 100
filepath = '../output/50_sq_degrees_100_patches/'
maps_count = 10

## TreeCorr
# for auto-correlation of the gaussian density fluctuations
# also multiply by mean density in patch to get integrated 3-pt function i_Xi

for j in range(maps_count):
    #print('gaussian Map # ',str(j+1))
    filepath_map = filepath+'gaussian_map_'+str(j+1)+'/'
    createFolder(filepath_map+'B_treecorr_patches_correlated/')

    for i in range(patch_count):
        #print('Patch # '+str(i+1))
        density_fluctuations = hp.read_cl(filepath_map+'A_healpy_patches_produced/del_gaussian_patch_'+str(i+1)+'.fits')

        density_fluctuations_RA = density_fluctuations[0,:]
        density_fluctuations_dec = density_fluctuations[1,:]
        density_fluctuations_del_g = density_fluctuations[2,:]
        
        mean_del_g = np.mean(density_fluctuations_del_g)
        var_del_g = np.var(density_fluctuations)

        cat = treecorr.Catalog(ra=density_fluctuations_RA, dec=density_fluctuations_dec,
                                ra_units='rad', dec_units='rad', k=density_fluctuations_del_g)

        kk = treecorr.KKCorrelation(min_sep=1, max_sep=400, nbins=30, sep_units='arcmin')
        kk.process(cat) 
        theta_tc = kk.meanr # tc stands for treecorr
        xi_tc = kk.xi
        mean_del_g_vec = np.ones(theta_tc.size)*mean_del_g
        var_del_g_vec = np.ones(theta_tc.size)*var_del_g
        i_Xi = xi_tc*mean_del_g # integrated 3-pt correlation function
        
        dat = np.array([theta_tc, xi_tc, mean_del_g_vec, var_del_g_vec, i_Xi]) # the mean_del_g_vec, var_del_g_vec is the same value for all angles (just for output purposes putting into column)
        dat = dat.T
        np.savetxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', dat, delimiter = ' ')


end = time.time()
print('\nTREECORR --> Time taken for execution (seconds): ', end - start)