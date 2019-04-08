import time
start = time.time()

import numpy as np
import healpy as hp
import os

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

# as l=0 and l=1 (and corresponding cl values of 0) are missing due to requirement of flask, we append them
# we also take only till l=8192 (before index 8191) for flask simulation
# therefore, we finally take l values from l=0 to l=8192 (and corresponding cl)

def read_cl():
    l = np.loadtxt('../../data/Cell_data-f1z1f1z1.dat', usecols=(0))
    l = np.append(np.array([0.0,1.0]), l[:8191])
    cl = np.loadtxt('../../data/Cell_data-f1z1f1z1.dat', usecols=(1))
    cl = np.append(np.array([0.0,0.0]), cl[:8191])
    return l, cl

l , cl = read_cl()

## HEALPY

maps_count = 10
NSIDE = 4096
createFolder('../simulations_output/gaussian_maps/')

for j in range(maps_count):    
    print('Gaussian Map # ',str(j+1))
    density_field_gaussian = hp.synfast(cl, NSIDE)
    hp.write_map("../simulations_output/gaussian_maps/gaussian_map_"+str(j+1)+".fits", density_field_gaussian, overwrite='true')

end = time.time()
print('\nTime taken for execution (seconds): ', end - start)
