import time, h5py, os, sys
import numpy as np
import matplotlib.pyplot as plt
import winsound

''' misc '''
def make_dir(path):
    if not os.path.exists(path):
        print('make path')
        os.mkdir(path)
    else:
        print('path exists')

def polyfit2d(x, y, z, kx=3, ky=3, order=None):
    # x, y: axis, 1d array
    # z: fitting data, 2d array
    # kx, ky: order to fit
    # order: kx + ky <= order, overwrites kx, ky statements

    xx, yy = np.meshgrid(x, y)
    coeffs = np.ones((kx+1, ky+1))

    a = np.zeros((coeffs.size, xx.size))

    # for each coefficient produce array x^i, y^j
    for index, (j, i) in enumerate(np.ndindex(coeffs.shape)):
        # do not include powers greater than order
        if order is not None and i + j > order:
            arr = np.zeros_like(xx)
        else:
            arr = coeffs[i, j] * xx**i * yy**j
        a[index] = arr.ravel()

    # do leastsq fitting and return leastsq result
    return np.linalg.lstsq(a.T, np.ravel(z), rcond=None)

def load_CXRO_data(E0):
	filename = os.path.join('../../LCLS/lcls_beamline_toolbox-beta/lcls_beamline_toolbox/xraybeamline2d/cxro_data/Be.csv')
	cxro_data = np.genfromtxt(filename, delimiter=',')
	energy = cxro_data[:,0]; delta = cxro_data[:,1]; beta = cxro_data[:,2]

	delta = np.interp(E0, energy, delta); beta = np.interp(E0, energy, beta)
	return delta, beta