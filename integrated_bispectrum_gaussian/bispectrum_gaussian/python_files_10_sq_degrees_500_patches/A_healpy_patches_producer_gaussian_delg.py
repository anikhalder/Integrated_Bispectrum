import time
start = time.time()

import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits
from astropy.table import Table, Column
import random
import os

def createFolder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

### Parameters to change according to patch size and count
# Make 500 patches (discs) of 10 sq. degree pixels
sq_degrees = 10
patch_radius = 0.031 #rad
patch_count = 500
filepath = '../output/10_sq_degrees_500_patches/'
maps_count = 10

hdul = fits.open('../output/gaussian_maps/gaussian_map_'+str(1)+'.fits')
NSIDE = hdul[1].header['NSIDE']
hdul.close

# Loop through the maps
for j in range(maps_count):
    #print('gaussian Map # ',str(j+1))
    density_field_gaussian = hp.fitsfunc.read_map('../output/gaussian_maps/gaussian_map_'+str(j+1)+'.fits')

    filepath_map = filepath+'gaussian_map_'+str(j+1)+'/'
    createFolder(filepath_map+'A_healpy_patches_produced/')
    createFolder(filepath_map+'C_plot_output/')

    plt.figure(figsize=(12,12))
    hp.mollview(density_field_gaussian, min=-1, max=4)
    hp.graticule()
    plt.savefig(filepath_map+'C_plot_output/density_field_gaussian_map_'+str(j+1)+'.pdf')

    ## HEALPY
    for i in range(patch_count):
        pixel_id = random.randint(1, density_field_gaussian.size - 1)
        #print('Patch # '+str(i+1)+': ', pixel_id)

        patch_center = hp.pix2vec(NSIDE, pixel_id)

        # find the pixels 
        pixels_indices_patch = hp.query_disc(NSIDE, patch_center, patch_radius)

        density_field_gaussian_patch = density_field_gaussian[pixels_indices_patch]
        DEC = np.pi/2 - hp.pix2ang(NSIDE, pixels_indices_patch)[0]
        RA = hp.pix2ang(NSIDE, pixels_indices_patch)[1]

        RA_data = Column(RA, name='RA', dtype='float') # column having ra
        DEC_data = Column(DEC, name='DEC', dtype='float') # column having DEC
        del_g_data = Column(density_field_gaussian_patch, name='del_g', dtype='float') # column having pixels_patch_values

        density_fluctuations_table = Table((RA_data, DEC_data, del_g_data))  
        density_fluctuations_table.write(filepath_map+'/A_healpy_patches_produced/del_gaussian_patch_'+str(i+1)+'.fits')

end = time.time()
print('\nHEALPY --> Time taken for execution (seconds): ', end - start)
