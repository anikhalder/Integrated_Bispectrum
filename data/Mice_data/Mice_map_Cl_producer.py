import numpy as np
import healpy as hp


number_field_mice = hp.fitsfunc.read_map('./mice_number_field.fits')
mean_number = np.mean(number_field_mice)
density_field_mice = (number_field_mice - mean_number) / mean_number

cl = hp.sphtfunc.anafast(density_field_mice)
l = np.arange(0,cl.size,1)

dat = np.array([l, cl])
dat = dat.T
np.savetxt('Cell_mice_map.dat', dat, delimiter = ' ')
