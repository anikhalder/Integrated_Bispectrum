import numpy as np

l = np.loadtxt('./Cell_CLASS.dat', usecols=(0,))
power_cl = np.loadtxt('./Cell_CLASS.dat', usecols=(1,))

cl = power_cl*2*np.pi/ l / (l+1)

dat = np.array([l, cl])
dat = dat.T
np.savetxt('Cell_data-f1z1f1z1.dat', dat, delimiter = ' ')
