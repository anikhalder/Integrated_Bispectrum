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

# as l=0 and l=1 (and corresponding cl values of 0) are missing due to requirement of flask, we append them
# we also take only till l=8192 (before index 8191) for flask simulation
# therefore, we finally take l values from l=0 to l=8192 (and corresponding cl)

def read_cl():
    l = np.loadtxt('../data/Cell_data-f1z1f1z1.dat', usecols=(0))
    l = np.append(np.array([0.0,1.0]), l[:8191])
    cl = np.loadtxt('../data/Cell_data-f1z1f1z1.dat', usecols=(1))
    cl = np.append(np.array([0.0,0.0]), cl[:8191])
    return l, cl

l , cl = read_cl()

# Theoretical angular correlation function (using formula with Legendre polynomials)
def w_theta(theta):
    x = np.cos(theta)
    coeff = (2*l+1)/(4*np.pi)*cl
    w = np.polynomial.legendre.legval(x, coeff)
    return w     

plt.figure(figsize=(9,9))

# -----------------------------------------------------------------------------------
# Compute mean of theta and i_Xi over all maps
theta_vec = np.loadtxt(filepath+'gaussian_map_1/B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(1)+'.txt', usecols=(0)) # in arcmins

theta_mean_all_maps_vec = np.zeros(theta_vec.size)
i_Xi_mean_all_maps_vec = np.zeros(theta_vec.size)

for j in range(maps_count):
    #print('gaussian Map # ',str(j+1))
    filepath_map = filepath+'gaussian_map_'+str(j+1)+'/'  

    theta_mean_one_map_vec = np.zeros(theta_vec.size)
    xi_mean_one_map_vec = np.zeros(theta_vec.size)
    mean_del_g_mean_one_map_vec = 0
    i_Xi_mean_one_map_vec = np.zeros(theta_vec.size)

    for i in range(patch_count):
        #print('Patch # '+str(i+1))
        # read individual map's treecorr data
        theta_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(0)) # in arcmins
        xi_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(1)) # xi aka w(theta)
        mean_del_g = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(2))[0]
        ##var_del_g = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(3))[0]
        i_Xi_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(4)) # i_Xi

        theta_mean_one_map_vec = theta_mean_one_map_vec + theta_vec
        xi_mean_one_map_vec = xi_mean_one_map_vec + xi_vec
        mean_del_g_mean_one_map_vec = mean_del_g_mean_one_map_vec + mean_del_g
        i_Xi_mean_one_map_vec = i_Xi_mean_one_map_vec + i_Xi_vec

    # individual map vector
    theta_mean_one_map_vec = theta_mean_one_map_vec/patch_count
    #i_Xi_mean_one_map_vec = i_Xi_mean_one_map_vec/patch_count 
    i_Xi_mean_one_map_vec = i_Xi_mean_one_map_vec/patch_count - mean_del_g_mean_one_map_vec/patch_count * xi_mean_one_map_vec/patch_count # for smaller errors

    # plot i_Xi of each map as a scatter plot
    plt.scatter(theta_mean_one_map_vec, i_Xi_mean_one_map_vec)

    theta_mean_all_maps_vec = theta_mean_all_maps_vec +theta_mean_one_map_vec
    i_Xi_mean_all_maps_vec = i_Xi_mean_all_maps_vec + i_Xi_mean_one_map_vec

# all maps mean vector
theta_mean_all_maps_vec = theta_mean_all_maps_vec/maps_count
i_Xi_mean_all_maps_vec = i_Xi_mean_all_maps_vec/maps_count

# -----------------------------------------------------------------------------------
# Compute variance of theta and i_Xi over all maps
theta_vec = np.loadtxt(filepath+'gaussian_map_1/B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(1)+'.txt', usecols=(0)) # in arcmins

theta_variance_all_maps_vec = np.zeros(theta_vec.size)
i_Xi_variance_all_maps_vec = np.zeros(theta_vec.size)

for j in range(maps_count):
    filepath_map = filepath+'gaussian_map_'+str(j+1)+'/'

    theta_mean_one_map_vec = np.zeros(theta_vec.size)
    xi_mean_one_map_vec = np.zeros(theta_vec.size)
    mean_del_g_mean_one_map_vec = 0
    i_Xi_mean_one_map_vec = np.zeros(theta_vec.size)

    for i in range(patch_count):
        #print('Patch # '+str(i+1))
        # read individual map's treecorr data
        theta_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(0)) # in arcmins
        xi_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(1)) # xi aka w(theta)
        mean_del_g = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(2))[0]
        ##var_del_g = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(3))[0]
        i_Xi_vec = np.loadtxt(filepath_map+'B_treecorr_patches_correlated/correlation_gaussian_patch_'+str(i+1)+'.txt', usecols=(4)) # i_Xi

        theta_mean_one_map_vec = theta_mean_one_map_vec + theta_vec
        xi_mean_one_map_vec = xi_mean_one_map_vec + xi_vec
        mean_del_g_mean_one_map_vec = mean_del_g_mean_one_map_vec + mean_del_g
        i_Xi_mean_one_map_vec = i_Xi_mean_one_map_vec + i_Xi_vec

    # individual map vector
    theta_mean_one_map_vec = theta_mean_one_map_vec/patch_count
    #i_Xi_mean_one_map_vec = i_Xi_mean_one_map_vec/patch_count
    i_Xi_mean_one_map_vec = i_Xi_mean_one_map_vec/patch_count - mean_del_g_mean_one_map_vec/patch_count * xi_mean_one_map_vec/patch_count # for smaller errors

    theta_variance_all_maps_vec = theta_variance_all_maps_vec + (theta_mean_one_map_vec - theta_mean_all_maps_vec)**2
    i_Xi_variance_all_maps_vec = i_Xi_variance_all_maps_vec + (i_Xi_mean_one_map_vec - i_Xi_mean_all_maps_vec)**2

# all maps variance vector
theta_variance_all_maps_vec = theta_variance_all_maps_vec/(maps_count-1)
i_Xi_variance_all_maps_vec = i_Xi_variance_all_maps_vec/(maps_count-1)

theta_std_dev_all_maps_vec = np.sqrt(theta_variance_all_maps_vec) # x-error
i_Xi_std_dev_all_maps_vec = np.sqrt(i_Xi_variance_all_maps_vec) # y-error

i_Xi_std_dev_mean_all_maps_vec = i_Xi_std_dev_all_maps_vec/np.sqrt(maps_count)

# -----------------------------------------------------------------------------------

createFolder(filepath+'plot_output/')

dat = np.array([theta_mean_all_maps_vec, i_Xi_mean_all_maps_vec, i_Xi_std_dev_all_maps_vec, i_Xi_std_dev_mean_all_maps_vec])
dat = dat.T
np.savetxt(filepath+'plot_output/i_Xi_simulations_gaussian_'+str(maps_count)+'_maps_'+str(patch_count)+'_patches_'+str(sq_degrees)+'_sq_degrees.txt', dat, delimiter = ' ')

#plt.scatter(theta_mean_one_map_vec, i_Xi_mean_vec, c='b', marker=10, label='bipectrum')
plt.errorbar(theta_mean_all_maps_vec, i_Xi_mean_all_maps_vec, yerr=i_Xi_std_dev_all_maps_vec, marker=10, label='i_Xi - one map error')
plt.errorbar(theta_mean_all_maps_vec, i_Xi_mean_all_maps_vec, yerr=i_Xi_std_dev_mean_all_maps_vec, marker=10, color='k', label='i_Xi - mean error')
#plt.plot(theta_mean_one_map_vec, w_theta(theta_mean_one_map_vec/60*np.pi/180), c='r', label='theoretical w(theta)')
plt.xlim(1,400)
#plt.ylim(1e-6, 1e-1)
plt.ylim(-0.001, 0.02)
plt.xscale('log')
#plt.yscale('log')
plt.axhline(0, linestyle='dashed')
plt.xlabel('Angle, theta (arcmins)', fontsize=14)
plt.ylabel('Integrated 3-pt function, i_Xi', fontsize=14)
plt.title('i_Xi of gaussian field - '+str(maps_count)+' maps with each map '+str(patch_count)+' patches ('+str(sq_degrees)+' sq degrees each patch)')
plt.legend(fontsize=13)
plt.savefig(filepath+'plot_output/i_Xi_gaussian_'+str(maps_count)+'_maps_'+str(patch_count)+'_patches_'+str(sq_degrees)+'_sq_degrees.pdf')

end = time.time()
print('\nANALYSIS --> Time taken for execution (seconds): ', end - start)