import time
start_program = time.time()

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import treecorr
import skinematics as skin
import os
import sys
import vegas

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
log_shift = 1.0 # log shift parameter
N = 10000 # number of evaluations for the integral
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
def w_cos_theta(cos_theta):
    x = cos_theta
    coeff = (2*l+1)/(4*np.pi)*cl
    w = np.polynomial.legendre.legval(x, coeff)
    return w     

def cos_angular_length(omega_1, omega_2):
    return np.cos(omega_1[0])*np.cos(omega_2[0]) + np.sin(omega_1[0])*np.sin(omega_2[0])*np.cos(omega_1[1]-omega_2[1])

def window_circular_patch(omega, omega_patch, patch_radius):
    if (np.arccos(cos_angular_length(omega, omega_patch)) <= patch_radius):
        return 1
    else:
        return 0

def bin_angular_scale_interval(omega_1, omega_2, theta_scale_interval):
    # theta_scale_interval is a list containing the lower and upper values of the bin (used in TreeCorr) in regular distance space in radians
    if (np.arccos(cos_angular_length(omega_1, omega_2)) >= theta_scale_interval[0] and np.arccos(cos_angular_length(omega_1, omega_2)) < theta_scale_interval[1]):
        return 1/(theta_scale_interval[1]-theta_scale_interval[0])
    else:
        return 0

theta_arr = np.linspace(0,2*patch_radius+0.001,10000)
cos_theta_arr = np.cos(theta_arr)
w_cos_theta_arr = w_cos_theta(cos_theta_arr)

w_cos_theta_interp = interpolate.interp1d(cos_theta_arr, w_cos_theta_arr)
    
def lognormal_3pt_corr(omega_1, omega_2, omega_3, log_shift):
    # angular 3pt correlartion function for a lognormal density field (according to Hilbert et al.)

    w_12 = w_cos_theta_interp(cos_angular_length(omega_1, omega_2))
    w_13 = w_cos_theta_interp(cos_angular_length(omega_1, omega_3))
    w_23 = w_cos_theta_interp(cos_angular_length(omega_2, omega_3))

    return (log_shift**(-1))*(w_12*w_13+w_12*w_23+w_13*w_23)+(log_shift**(-3))*(w_12*w_13*w_23)

"""
#########################
### Scipy integration ###
#########################

#### DOES NOT CONVERGE ####

def calculate_patch_area_scipy(patch_radius):
    # patch is basically a spherical cap centered around the North Pole one a sphere of unit radius
    def integrand(theta, phi):
        return np.sin(theta)

    area = integrate.nquad(integrand, [[0, patch_radius], [0, 2*np.pi]])[0]

    return area

def integrated_lognormal_3pt_corr_scipy(theta_scale_interval, patch_radius, log_shift, patch_area):

    def integrand(theta_1, phi_1, theta_2, phi_2, theta_3, phi_3):
        return np.sin(theta_1)*np.sin(theta_2)*np.sin(theta_3)*lognormal_3pt_corr([theta_1, phi_1], [theta_2, phi_2], [theta_3, phi_3], log_shift)*bin_angular_scale_interval([theta_1, phi_1], [theta_3, phi_3], theta_scale_interval)

    i_Xi = integrate.nquad(integrand, [[0, patch_radius], [0, 2*np.pi], [0, patch_radius], [0, 2*np.pi], [0, patch_radius], [0, 2*np.pi]])

    return i_Xi[0] / (patch_area**2)
"""

#########################
### Vegas integration ###
#########################

def calculate_patch_area_vegas(patch_radius):
    # patch is basically a spherical cap centered around the North Pole one a sphere of unit radius
    def integrand(x): # x[0] stands for theta, x[1] stands for phi
        return np.sin(x[0]) 

    integ = vegas.Integrator([[0, patch_radius], [0, 2*np.pi]])
    area = integ(integrand, nitn=20, neval=5000)
    #print(result.summary())
    #print('result = %s    Q = %.2f' % (result, result.Q))
    
    return area.mean  # return the weighted mean over the total number of iterations

def integrated_lognormal_3pt_corr_vegas(theta_scale_interval, patch_radius, log_shift, patch_area, theta_scale, N=10000):

    def integrand(x):
        theta_1 = x[0]
        phi_1 = x[1]
        theta_2 = x[2]
        phi_2 = x[3]
        theta_3 = x[4]
        phi_3 = x[5]
    
        return np.sin(theta_1)*np.sin(theta_2)*lognormal_3pt_corr([theta_1, phi_1], [theta_2, phi_2], [theta_3, phi_3], log_shift)*bin_angular_scale_interval([theta_1, phi_1], [theta_3, phi_3], theta_scale_interval)

    integrator = vegas.Integrator([[0, patch_radius], [0, 2*np.pi], [0, patch_radius], [0, 2*np.pi], [0, patch_radius], [0, 2*np.pi]])

    i_Xi = integrator(integrand, nitn=20, neval=N)

    return i_Xi.mean / (2*np.pi*np.sin(theta_scale)) / (patch_area**2) # return the weighted mean over the total number of iterations divided by the patch_area squared


##########################
### Manual integration ###
##########################

def spherical_to_cartesian(P_sph, sph_units):
    # r , theta, phi 
    r = P_sph[0]
    theta = 0
    phi = 0
    if (sph_units == 'deg'):
        theta = np.deg2rad(P_sph[1])
        phi = np.deg2rad(P_sph[2])
    elif (sph_units == 'rad'):
        theta = P_sph[1]
        phi = P_sph[2]  
    else:
        raise Exception('units should be either \'deg\' or \'rad\'')
   
    return np.array([r*np.sin(theta)*np.cos(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(theta)])

def cartesian_to_spherical(P_cart, sph_units):
    # x, y, z ---> r, theta, phi 
    x = P_cart[0]
    y = P_cart[1]
    z = P_cart[2]
    
    r = np.sqrt(x*x+y*y+z*z)
    
    theta = np.arctan2(np.sqrt(x*x+y*y), z)
        
    phi = np.arctan2(y,x)
    if phi < 0:
        phi = 2*np.pi + phi
    
    if (sph_units == 'deg'):
        return np.array([r, theta*180/np.pi, phi*180/np.pi])
    elif (sph_units == 'rad'):
        return np.array([r, theta, phi])
    else:
        raise Exception('units should be either \'deg\' or \'rad\'')

def rotate_vector_to_NP_matrices(P_sph, P_sph_units):
    # retrieve matrices to rotate a vector to the North Pole
    theta = 0
    phi = 0
    if (P_sph_units == 'deg'):
        theta = P_sph[1]
        phi = P_sph[2]         
    elif (P_sph_units == 'rad'):
        theta = np.rad2deg(P_sph[1])
        phi = np.rad2deg(P_sph[2])
    else:
        raise Exception('units should be either \'deg\' or \'rad\'')
        
    # First rotate about z axis (counter azimuthally) so as to make the point lie in the xz plane i.e. y = 0
    Rot_z = skin.rotmat.R(axis='z', angle=-phi) # skin.rotmar.R only accepts angle in degrees!   
    
    # Then rotate about y axis (counter polar-wise) to make the point lie at the NP i.e. x = 0, y = 0
    Rot_y = skin.rotmat.R(axis='y', angle=-theta) # skin.rotmar.R only accepts angle in degrees!
       
    return np.array([Rot_z, Rot_y]) # 1st rotation matrix and then 2nd rotation matrix

def draw_pt_within_patch(patch_radius, patch_area):

    theta = np.random.uniform(0, patch_radius, 1)[0]
    phi = np.random.uniform(0, 2*np.pi, 1)[0]
    
    # random point in the patch using the ideology behind picking a point uniformly on a unit sphere from http://corysimon.github.io/articles/uniformdistn-on-sphere/
    #theta = np.arccos(1-patch_area*np.random.uniform(0, 1, 1)[0]/(2*np.pi)) # Anik transformation
    #theta = np.arccos(np.cos(patch_radius)+patch_area*np.random.uniform(0, 1, 1)[0]/(2*np.pi)) # Oliver transformation

    return [theta, phi]

def draw_pt_on_circle_within_patch(omega_1, theta_scale, patch_radius):
    R_1, R_2 = rotate_vector_to_NP_matrices([1, omega_1[0], omega_1[1]], 'rad')
    
    theta_2_prime = theta_scale # polar angle in the rotated frame
    
    omega_2_3d_sph = np.zeros(3)
    phi_2_prime = np.random.uniform(0, 2*np.pi, 1)[0]
    omega_2_prime_cart = spherical_to_cartesian([1, theta_2_prime, phi_2_prime], 'rad') # set radius of sphere to 1
    omega_2_3d_sph = cartesian_to_spherical(np.matmul(R_1.T, np.matmul(R_2.T, omega_2_prime_cart)), 'rad')

    if (omega_2_3d_sph[1] > patch_radius):
        return [float('nan'), omega_2_3d_sph[2]]
    else:
        return [omega_2_3d_sph[1], omega_2_3d_sph[2]]

def integrated_lognormal_3pt_corr_manual(theta_scale, patch_radius, log_shift, patch_area, N=10000):
    f = 0
    accepted_points = 0
    for i in range(N):
        theta_1, phi_1 = draw_pt_within_patch(patch_radius, patch_area)
        theta_2, phi_2 = draw_pt_within_patch(patch_radius, patch_area)
        theta_3, phi_3 = draw_pt_on_circle_within_patch([theta_1, phi_1], theta_scale, patch_radius)

        if (not np.isnan(theta_3)):
            accepted_points += 1
            # using uniform theta and uniform phi in draw_pt_within_patch()       
            f += np.sin(theta_1)*np.sin(theta_2)*lognormal_3pt_corr([theta_1, phi_1], [theta_2, phi_2], [theta_3, phi_3], log_shift)

            # using transformed theta and uniform phi in draw_pt_within_patch()          
            #f += lognormal_3pt_corr([theta_1, phi_1], [theta_2, phi_2], [theta_3, phi_3], log_shift)

    # using uniform theta and uniform phi in draw_pt_within_patch()     
    return f / accepted_points * (2*np.pi*patch_radius/patch_area)**2

    # using transformed theta and uniform phi in draw_pt_within_patch()
    #return f / accepted_points

#######################################################################################################################################
#######################################################################################################################################
# Individual function tests

#t1 = [0.031, 3.2]
t1 = [0.003, 0]
t2 = [2.1, 0.7]
t3 = [0.0155, 3.141592653589793]
t4 = [0.010936910628127367, 3.141592653589793]
t5 = [0.020063089371872633, 3.141592653589793]
t6 = [0.0155, 3.141592653589793]
t7 = [0.0312, 0.4]
t8 = [0.0309, 0.4]
tL = [0.0, 0.0] # location of patch -> at the north pole

t_scale = 0.003 # test scale we are interested (in radians)

print("Angular length between t1 and t2 (in radians) = ", cos_angular_length(t1, t2))
print("Angular length between t1 and tL (in radians)  = ", cos_angular_length(t1, tL))
print("Angular length between t2 and tL (in radians)  = ", cos_angular_length(t2, tL))
print("Angular length between t3 and tL (in radians)  = ", cos_angular_length(t3, tL))
print("Angular length between t4 and t6 (in radians)  = ", cos_angular_length(t4, t6))
print("Angular length between t5 and t6 (in radians)  = ", cos_angular_length(t5, t6))

print("Is t1 within the patch which is centered at tL and is of size "+str(patch_radius)+" radians = ", window_circular_patch(t1, tL, patch_radius))
print("Is t2 within the patch which is centered at tL and is of size "+str(patch_radius)+" radians = ", window_circular_patch(t2, tL, patch_radius))
print("Is t7 within the patch which is centered at tL and is of size "+str(patch_radius)+" radians = ", window_circular_patch(t7, tL, patch_radius))
print("Is t8 within the patch which is centered at tL and is of size "+str(patch_radius)+" radians = ", window_circular_patch(t8, tL, patch_radius))

print("Is t1 and tL at approximately "+str(t_scale)+" radians separation =", bin_angular_scale_interval(t1, tL, [t_scale - 0.00001, t_scale + 0.00001]))

print("Angular 2-pt correlation: w(theta between t1 and t2) = ", w_cos_theta(cos_angular_length(t1, t2)))
print("Angular 2-pt correlation: w(theta between t4 and t6) = ", w_cos_theta(cos_angular_length(t4, t6)))
print("Angular 2-pt correlation: w(theta between t5 and t6) = ", w_cos_theta(cos_angular_length(t5, t6)))

#print("Lognormal 3-pt correlation between t1, t2, t3 = ", lognormal_3pt_corr(t1, t2, t3, log_shift))

"""
print("\n #### Scipy integration ###")

A_L = calculate_patch_area(patch_radius)
print("Area of the patch = ", A_L)

print("Integrated lognormal 3-pt function at angular scale "+str(t_scale)+" radians = ", integrated_lognormal_3pt_corr_scipy([t_scale - 0.00001, t_scale + 0.00001], t_patch_radius, log_shift, A_L))
"""


print("\n #### Vegas integration ###")

A_L = calculate_patch_area(patch_radius)
print("Area of the patch = ", A_L)

#print("Integrated lognormal 3-pt function at angular scale "+str(t_scale)+" radians = ", integrated_lognormal_3pt_corr_vegas([t_scale - 0.00001, t_scale + 0.00001], t_patch_radius, log_shift, A_L))

"""
print("\n #### Manual integration ###")

A_L = calculate_patch_area(patch_radius)
print("Area of the patch = ", A_L)
print("Integrated lognormal 3-pt function at angular scale "+str(t_scale)+" radians = ", integrated_lognormal_3pt_corr_manual(t_scale, patch_radius, log_shift, A_L))
"""

#######################################################################################################################################
#######################################################################################################################################
# Evaluate the integrated bispectrum (for lognormal density field) for given scale angles (for Vegas - with same bin intervals as used by treecorr for doing the numerical work with the patches)

kk = treecorr.KKCorrelation(min_sep=1, max_sep=400, nbins=30, sep_units='arcmin')
theta_scale_log_vec = kk.logr # in log arcmins

#########################
### Vegas integration ###
#########################

createFolder(str(sq_degrees)+'_sq_degrees/'+str(N)+'_pts/')

bin_width_log = theta_scale_log_vec[1]-theta_scale_log_vec[0]

theta_scale_log_bin_edges = np.zeros(theta_scale_log_vec.size + 1)
theta_scale_log_bin_edges[0] = theta_scale_log_bin_edges[0]-bin_width_log/2# in log arcmins
for j in range(len(theta_scale_log_vec)):
    theta_scale_log_bin_edges[j+1] = theta_scale_log_vec[j]+bin_width_log/2 # in log arcmins

theta_scale_bin_edges = np.exp(theta_scale_log_bin_edges) # in arcmins

theta_scale_vec = np.exp(theta_scale_log_vec) # in arcmins
i_Xi_vec = np.zeros(theta_scale_log_vec.size)
for i in range(theta_scale_log_vec.size):
    print('###############################################################################')
    print('Iteration #',str(i+1))
    start_iter = time.time()
    print('Theta (in arcmins): ', theta_scale_vec[i]) # nominal center of the bin used by treecorr
    theta_scale_interval = [theta_scale_bin_edges[i]*np.pi/180/60, theta_scale_bin_edges[i+1]*np.pi/180/60]  # from arcmins to radians
    i_Xi_vec[i] = integrated_lognormal_3pt_corr_vegas(theta_scale_interval, patch_radius, log_shift, A_L, theta_scale_vec[i], N)
    print('i_Xi: ', i_Xi_vec[i])
    end_iter = time.time()
    print('Time taken for execution of iteration (seconds): ', end_iter - start_iter)

"""
##########################
### Manual integration ###
##########################

createFolder(str(sq_degrees)+'_sq_degrees/'+str(N)+'_pts/')

theta_scale_vec = np.exp(theta_scale_log_vec)
i_Xi_vec = np.zeros(theta_scale_vec.size)
for i in range(theta_scale_vec.size):
    # the last 3/5 theta scales are bigger than the patch radius for 10/5 sq degrees
    print('###############################################################################')
    print('Iteration #',str(i+1))
    start_iter = time.time()
    print('Theta (in arcmins): ', theta_scale_vec[i]) # nominal center of the bin used by treecorr
    i_Xi_vec[i] = integrated_lognormal_3pt_corr_manual(theta_scale_vec[i]*np.pi/180/60, patch_radius, log_shift, A_L, N)
    print('i_Xi: ', i_Xi_vec[i])
    end_iter = time.time()
    print('Time taken for execution of iteration (seconds): ', end_iter - start_iter)
"""

dat = np.array([theta_scale_vec, i_Xi_vec])

dat = dat.T
np.savetxt(str(sq_degrees)+'_sq_degrees/'+str(N)+'_pts/i_Xi_theoretical_lognormal_patch_'+str(sq_degrees)+'_sq_degrees.txt', dat, delimiter = ' ')    

####################################
# Theoretical i_Xi plot
####################################
theta_scale_vec = np.loadtxt(str(sq_degrees)+'_sq_degrees/'+str(N)+'_pts/i_Xi_theoretical_lognormal_patch_'+str(sq_degrees)+'_sq_degrees.txt', usecols=(0)) # in arcmins
i_Xi_vec = np.loadtxt(str(sq_degrees)+'_sq_degrees/'+str(N)+'_pts/i_Xi_theoretical_lognormal_patch_'+str(sq_degrees)+'_sq_degrees.txt', usecols=(1))
plt.figure(figsize=(9,9))
plt.plot(theta_scale_vec, i_Xi_vec, c='r', label='theoretical i_Xi(theta)')
plt.xlim(1,400)
#plt.ylim(1e-6, 1e-1)
plt.ylim(-0.001, 0.02)
plt.xscale('log')
#plt.yscale('log')
plt.axhline(0, linestyle='dashed')
plt.xlabel('Angle, theta (arcmins)', fontsize=14)
plt.ylabel('Integrated 3-pt function, i_Xi', fontsize=14)
plt.title('Integrated 3-pt function of lognormal field ('+str(sq_degrees)+' sq degrees patch) \n # of evaluations = '+str(N))
plt.legend(fontsize=13)
plt.savefig(str(sq_degrees)+'_sq_degrees/'+str(N)+'_pts/i_Xi_theoretical_lognormal_patch_'+str(sq_degrees)+'_sq_degrees.pdf')

####################################
# Theoretical & simulations i_Xi plot
####################################
theta_mean_all_maps_vec = np.loadtxt(str(sq_degrees)+'_sq_degrees/i_Xi_simulations_lognormal_'+str(maps_count)+'_maps_'+str(patch_count)+'_patches_'+str(sq_degrees)+'_sq_degrees.txt', usecols=(0)) # in arcmins
i_Xi_mean_all_maps_vec = np.loadtxt(str(sq_degrees)+'_sq_degrees/i_Xi_simulations_lognormal_'+str(maps_count)+'_maps_'+str(patch_count)+'_patches_'+str(sq_degrees)+'_sq_degrees.txt', usecols=(1))
i_Xi_std_dev_all_maps_vec = np.loadtxt(str(sq_degrees)+'_sq_degrees/i_Xi_simulations_lognormal_'+str(maps_count)+'_maps_'+str(patch_count)+'_patches_'+str(sq_degrees)+'_sq_degrees.txt', usecols=(2))
i_Xi_std_dev_mean_all_maps_vec = np.loadtxt(str(sq_degrees)+'_sq_degrees/i_Xi_simulations_lognormal_'+str(maps_count)+'_maps_'+str(patch_count)+'_patches_'+str(sq_degrees)+'_sq_degrees.txt', usecols=(3))


plt.figure(figsize=(9,9))
plt.plot(theta_scale_vec, i_Xi_vec, c='r', label='theoretical i_Xi(theta)')
plt.errorbar(theta_mean_all_maps_vec, i_Xi_mean_all_maps_vec, yerr=i_Xi_std_dev_all_maps_vec, marker=10, label='i_Xi - one map error')
plt.errorbar(theta_mean_all_maps_vec, i_Xi_mean_all_maps_vec, yerr=i_Xi_std_dev_mean_all_maps_vec, marker=10, color='k', label='i_Xi - mean error')
plt.xlim(1,400)
#plt.ylim(1e-6, 1e-1)
plt.ylim(-0.001, 0.02)
plt.xscale('log')
#plt.yscale('log')
plt.axhline(0, linestyle='dashed')
plt.xlabel('Angle, theta (arcmins)', fontsize=14)
plt.ylabel('Integrated 3-pt function, i_Xi', fontsize=14)
plt.title('Integrated 3-pt function of lognormal field ('+str(sq_degrees)+' sq degrees patch) \n # of evaluations = '+str(N))
plt.legend(fontsize=13)
plt.savefig(str(sq_degrees)+'_sq_degrees/'+str(N)+'_pts/i_Xi_lognormal_theoretical_simulations_patch_'+str(sq_degrees)+'_sq_degrees.pdf')

####################################
# Ratio plot
####################################
plt.figure(figsize=(9,9))
plt.scatter(theta_mean_all_maps_vec, i_Xi_mean_all_maps_vec/i_Xi_vec)
plt.xlim(1,400)
plt.ylim(0,5)
plt.xscale('log')
plt.axhline(1, linestyle='dashed')
plt.xlabel('Angle, theta (arcmins)', fontsize=14)
plt.ylabel('Ratio', fontsize=14)
plt.title('Ratio of simulation and theoretical i_Xi ('+str(sq_degrees)+' sq degrees patch) \n # of evaluations = '+str(N))
plt.savefig(str(sq_degrees)+'_sq_degrees/'+str(N)+'_pts/i_Xi_lognormal_ratio_simulations_theoretical_patch_'+str(sq_degrees)+'_sq_degrees.pdf')


end_program = time.time()
print('\nAnalysis --> Time taken for execution (seconds): ', end_program - start_program)