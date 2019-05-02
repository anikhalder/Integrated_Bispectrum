import time
start = time.time()

import numpy as np
import matplotlib.pyplot as plt
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
    #filepath = '../mice_output/250_sq_degrees_20_patches/'

elif (sys.argv[1] == str(50)):
    ### Parameters to change according to patch size and count
    # Make 100 patches (discs) of 50 sq. degree pixels
    sq_degrees = 50
    #patch_radius = 0.069 #rad
    patch_count = 100
    #filepath = '../mice_output/50_sq_degrees_100_patches/'

elif (sys.argv[1] == str(10)):
    ### Parameters to change according to patch size and count
    # Make 500 patches (discs) of 10 sq. degree pixels
    sq_degrees = 10
    #patch_radius = 0.031 #rad
    patch_count = 500
    #filepath = '../mice_output/10_sq_degrees_500_patches/'

elif (sys.argv[1] == str(5)):
    ### Parameters to change according to patch size and count
    # Make 1000 patches (discs) of 5 sq. degree pixels
    sq_degrees = 5
    #patch_radius = 0.022 #rad
    patch_count = 1000
    #filepath = '../mice_output/5_sq_degrees_1000_patches/'
    
else:
	raise Exception('Choose correct patch size!')

def calculate_patch_radius(patch_area_sq_degrees):
    return np.arccos(1-patch_area_sq_degrees*np.pi/(2*180*180))

def calculate_patch_area(patch_radius): 
    # patch is basically a spherical cap centered around the North Pole on a sphere of unit radius
    r = 1
    return 2*np.pi*r*r*(1-np.cos(patch_radius)) # from Wikipedia

patch_radius = calculate_patch_radius(float(sys.argv[1]))
filepath = '../mice_output/'+str(sq_degrees)+'_sq_degrees_'+str(patch_count)+'_patches/'

################################################################

plt.figure(figsize=(10,10))

# Compute mean of theta and i_zeta over all patches in the mice map
theta_vec = np.loadtxt(filepath+'mice_map/B_treecorr_patches_correlated/correlation_mice_patch_'+str(1)+'.txt', usecols=(0)) # in arcmins

filepath_map = filepath+'mice_map/'  

theta_mean_mice_map_vec = np.zeros(theta_vec.size)
w_mean_mice_map_vec = np.zeros(theta_vec.size)
mean_del_mean_mice_map_vec = 0
i_zeta_mean_mice_map_vec = np.zeros(theta_vec.size)

for i in range(patch_count):
    # read mice map's treecorr data
    theta_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_mice_patch_'+str(i+1)+'.txt', usecols=(0)) # in arcmins
    w_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_mice_patch_'+str(i+1)+'.txt', usecols=(1)) # w(theta)
    mean_del = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_mice_patch_'+str(i+1)+'.txt', usecols=(2))[0]
    ##var_del = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_mice_patch_'+str(i+1)+'.txt', usecols=(3))[0]
    i_zeta_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_mice_patch_'+str(i+1)+'.txt', usecols=(4)) # i_zeta

    theta_mean_mice_map_vec = theta_mean_mice_map_vec + theta_vec
    w_mean_mice_map_vec = w_mean_mice_map_vec + w_vec
    mean_del_mean_mice_map_vec = mean_del_mean_mice_map_vec + mean_del
    i_zeta_mean_mice_map_vec = i_zeta_mean_mice_map_vec + i_zeta_vec

# individual map vector
theta_mean_mice_map_vec = theta_mean_mice_map_vec/patch_count
#i_zeta_mean_mice_map_vec = i_zeta_mean_mice_map_vec/patch_count 
i_zeta_mean_mice_map_vec = i_zeta_mean_mice_map_vec/patch_count - mean_del_mean_mice_map_vec/patch_count * w_mean_mice_map_vec/patch_count # for smaller errors 

# plot i_zeta of mice map as a scatter plot
#plt.scatter(theta_mean_mice_map_vec[1:], i_zeta_mean_mice_map_vec[1:])

# -----------------------------------------------------------------------------------
# Compute variance of i_zeta over all patches in the mice map
theta_vec = np.loadtxt(filepath+'mice_map/B_treecorr_patches_correlated/correlation_mice_patch_'+str(1)+'.txt', usecols=(0)) # in arcmins

i_zeta_variance_mice_map_vec = np.zeros(theta_vec.size)

filepath_map = filepath+'mice_map/'

for i in range(patch_count):
    i_zeta_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_mice_patch_'+str(i+1)+'.txt', usecols=(4)) # i_zeta

    i_zeta_variance_mice_map_vec = i_zeta_variance_mice_map_vec + (i_zeta_mean_mice_map_vec - i_zeta_vec)**2

i_zeta_variance_mice_map_vec = i_zeta_variance_mice_map_vec/(patch_count-1)
i_zeta_std_dev_mice_map_vec = np.sqrt(i_zeta_variance_mice_map_vec) # y-error
i_zeta_std_dev_mean_mice_map_vec = i_zeta_std_dev_mice_map_vec / np.sqrt(patch_count) # y-error

# -----------------------------------------------------------------------------------

createFolder(filepath+'plot_output/')

dat = np.array([theta_mean_mice_map_vec, i_zeta_mean_mice_map_vec, i_zeta_std_dev_mice_map_vec, i_zeta_std_dev_mean_mice_map_vec])
dat = dat.T
np.savetxt(filepath+'plot_output/i_zeta_mice_map_'+str(patch_count)+'_patches_'+str(sq_degrees)+'_sq_degrees.txt', dat, delimiter = ' ')

plt.errorbar(theta_mean_mice_map_vec[1:], i_zeta_mean_mice_map_vec[1:], yerr=i_zeta_std_dev_mice_map_vec[1:], marker=10, label='$i\\zeta$ - one patch error')
plt.errorbar(theta_mean_mice_map_vec[1:], i_zeta_mean_mice_map_vec[1:], yerr=i_zeta_std_dev_mean_mice_map_vec[1:], marker=10,  color='k', label='$i\\zeta$ - mean patches error')
plt.xlim(1,200)
#plt.ylim(1e-6, 1e-1)
plt.ylim(-0.001, 0.02)
plt.xscale('log')
#plt.yscale('log')
plt.axhline(0, linestyle='dashed')
plt.xlabel('Angle, $\\theta$ (arcmins)', fontsize=16)
plt.ylabel('Integrated 3-pt function, $i\\zeta(\\theta)$', fontsize=16)
plt.tick_params(labelsize=16)
plt.title('$i\\zeta$ of mice map - '+str(patch_count)+' patches ('+str(sq_degrees)+' sq. degrees each patch)', fontsize=14)
plt.legend(loc='best', fontsize=14)
plt.savefig(filepath+'plot_output/i_zeta_mice_map_'+str(patch_count)+'_patches_'+str(sq_degrees)+'_sq_degrees.pdf')

end = time.time()
print('\nANALYSIS --> Time taken for execution (seconds): ', end - start)