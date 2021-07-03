import time, h5py, os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from lcls_beamline_toolbox.xraybeamline2d import beam1d as beam, optics1d as optics, beamline1d as beamline
from lcls_beamline_toolbox.xraybeamline2d.util import Util

''' misc '''
def make_dir(path):
    if not os.path.exists(path):
        print('make path')
        os.mkdir(path)
    else:
        print('path exists')

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

def load_crystal_data(num, dir_profile):
    coords = np.loadtxt(dir_profile+'Nlist_%d.txt' %num, skiprows=12)
    data = np.loadtxt(dir_profile+'Uy_list_%d.txt' % num, skiprows=0)

    x = coords[:,1]
    y = coords[:,2]
    z = coords[:,3]
    dy = data[:,1]
    print(data.shape)
    # plt.figure()
    # plt.plot(x)
    # plt.figure()
    # plt.plot(z)

    xx = np.linspace(np.min(x), np.max(x), 1024)
    zz = np.linspace(np.min(z), np.max(z), 1024)
    xx,zz = np.meshgrid(xx,zz)

    # f = interpolate.interp2d(z,x,dy)
    dy2 = interpolate.griddata((x,z),dy,(xx,zz), fill_value=0)

    dy_symmetrize = np.concatenate((dy2,np.flipud(dy2)),axis=0)
    dy_symmetrize = np.concatenate((dy_symmetrize,np.fliplr(dy_symmetrize)),axis=1)

    xx2 = np.linspace(np.min(x), -np.min(x), 2048)
    zz2 = np.linspace(np.min(z), -np.min(z), 2048)
    xx2,zz2 = np.meshgrid(xx2,zz2)
    
    return dy_symmetrize, xx2, zz2

def find_zero(x, y):
    # if linearly decrease:
    if y[:5].mean() >= y[-5:].mean():
        index = np.where(y>=0)[0][-1]
    else:
        index = np.where(y<=0)[0][-1]
    slope = (y[index+1] - y[index])/(x[index+1] - x[index])
    result = x[index] - y[index]/slope
    print('left {}, right {}, result {}'.format(x[index], x[index+1], result))
    return result

''' define beamline '''
def define_Telescope(E0, m1_p=185.0, m2_p=175.5):
    z_s = 650
    
    ## Telescope
    m1 = optics.CurvedMirror('M1', p=m1_p, q=-58, length=1, z=185+z_s, alpha=2.65e-3)
    im_after_T1 = optics.PPM('im_after_T1', z=m1.z+.01, FOV=5e-3, N=256)
    
    m2 = optics.CurvedMirror('M2', p=m2_p, q=1e5, length=1, z=300+z_s, alpha=2.65e-3, orientation=2)
    im_after_T2 = optics.PPM('im_after_T2', z=m2.z+.01, FOV=5e-3, N=256)

    Telescope_devices = [m1, im_after_T1, m2, im_after_T2]

    return Telescope_devices

def define_HHLM_2DCM(
    E0, HHLM_distance_factor=1.0,
    hkl1 = [1,1,1], alphaAsym1 = 9.0,
    hkl2 = [1,1,1], alphaAsym2 = 0.0,
    shapeErrors=[None for i in range(6)],
    l_crystal=[1e-1 for i in range(6)],
    w_crystal = [5e-3 for i in range(6)]):
    
    '''
    defines the HHLM optics for the 2DCM setup (1-1-2-2)
    E0: photon energy [eV]
    hkl: crystal reflection surface indices for pair 1 and pair 2
    alphaAsym: asymmetry angle [degree], negative value means beam footprint increases after reflection
    shapeErrors: crytal shapeError as loaded from Lin's profiles

    returns optics
    '''
    z_s = 650
    
    ## HHLM
    asym1 = np.deg2rad(alphaAsym1)
    asym2 = np.deg2rad(alphaAsym2)
    hhlm1 = optics.Crystal('HHLM1', hkl=hkl1, length=l_crystal[0], width=w_crystal[0],
                           z=305+z_s, alphaAsym=-asym1, E0=E0, orientation=0, pol='s',
                           shapeError=shapeErrors[0])
    im_after_HHLM1 = optics.PPM('im_after_HHLM1', FOV=20e-3,N=256,z=hhlm1.z+.01)

    hhlm2 = optics.Crystal('HHLM2', hkl=hkl1, length=l_crystal[3], width=w_crystal[3],
                           z=hhlm1.z+139e-3*HHLM_distance_factor, alphaAsym=asym1, E0=E0, orientation=2, pol='s',
                           shapeError=shapeErrors[3])
    im_after_HHLM2 = optics.PPM('im_after_HHLM2', FOV=20e-3,N=256,z=hhlm2.z+.01)

    hhlm3 = optics.Crystal('HHLM3', hkl=hkl2, length=l_crystal[1], width=w_crystal[1],
                           z=hhlm1.z+361e-3*HHLM_distance_factor, alphaAsym=-asym2, E0=E0, orientation=2, pol='s',
                           shapeError=shapeErrors[1])
    im_after_HHLM3 = optics.PPM('im_after_HHLM3', FOV=20e-3,N=256,z=hhlm3.z+.01)

    hhlm4 = optics.Crystal('HHLM4', hkl=hkl2, length=l_crystal[2], width=w_crystal[2],
                           z=hhlm1.z+.5*HHLM_distance_factor, alphaAsym=asym2, E0=E0, orientation=0, pol='s',
                           shapeError=shapeErrors[2])
    im_after_HHLM4 = optics.PPM('im_after_HHLM4', FOV=20e-3,N=256,z=hhlm4.z+.01)

    HHLM_devices = [hhlm1, im_after_HHLM1, hhlm2, im_after_HHLM2, hhlm3, im_after_HHLM3, hhlm4, im_after_HHLM4]

    return HHLM_devices

def define_HHLM_Zigzag(
    E0, HHLM_distance_factor=1.0,
    hkl1 = [1,1,1], alphaAsym1 = 9.0,
    hkl2 = [4,4,0], alphaAsym2 = 0.0,
    shapeErrors=[None for i in range(6)],
    l_crystal=[1e-1 for i in range(6)],
    w_crystal = [5e-3 for i in range(6)]):
    '''
    defines the HHLM optics for the zigzag setup (1-2-2-1)
    E0: photon energy [eV]
    hkl: crystal reflection surface indices for pair 1 and pair 2
    alphaAsym: asymmetry angle [degree], negative value means beam footprint increases after reflection
    shapeErrors: crytal shapeError as loaded from Lin's profiles

    returns optics
    '''
    z_s = 650

    ## HHLM
    asym1 = np.deg2rad(alphaAsym1)
    asym2 = np.deg2rad(alphaAsym2)
    hhlm1 = optics.Crystal('HHLM1', hkl=hkl1, length=l_crystal[0], width=w_crystal[0],
                           z=305+z_s, alphaAsym=-asym1, E0=E0, orientation=0, pol='s',
                           shapeError=shapeErrors[0])
    im_after_HHLM1 = optics.PPM('im_after_HHLM1', FOV=30e-3,N=256,z=hhlm1.z+.01)

    hhlm2 = optics.Crystal('HHLM2', hkl=hkl2, length=l_crystal[1], width=w_crystal[1],
                           z=hhlm1.z+139e-3*HHLM_distance_factor, alphaAsym=-asym2, E0=E0,orientation=2, pol='s',
                           shapeError=shapeErrors[1])
    im_after_HHLM2 = optics.PPM('im_after_HHLM2', FOV=30e-3,N=256,z=hhlm2.z+.01)

    hhlm3 = optics.Crystal('HHLM3', hkl=hkl2, length=l_crystal[2], width=w_crystal[2],
                           z=hhlm1.z+361e-3*HHLM_distance_factor, alphaAsym=asym2, E0=E0,orientation=0, pol='s',
                           shapeError=shapeErrors[2])
    im_after_HHLM3 = optics.PPM('im_after_HHLM3', FOV=30e-3,N=256,z=hhlm3.z+.01)

    hhlm4 = optics.Crystal('HHLM4', hkl=hkl1, length=l_crystal[3], width=w_crystal[3],
                           z=hhlm1.z+.5*HHLM_distance_factor, alphaAsym=asym1, E0=E0,orientation=2, pol='s',
                           shapeError=shapeErrors[3])
    im_after_HHLM4 = optics.PPM('im_after_HHLM4', FOV=20e-3,N=256,z=hhlm4.z+.01)

    HHLM_devices = [hhlm1, im_after_HHLM1, hhlm2, im_after_HHLM2, hhlm3, im_after_HHLM3, hhlm4, im_after_HHLM4]

    return HHLM_devices

def define_HRM(E0, f1=10., f2=10., slit_width=3e-6,
               hkl=[4,4,0],alphaAsym=15.0,
               shapeErrors=[None for i in range(6)],
               l_crystal=[1e-1 for i in range(6)],
               w_crystal = [5e-3 for i in range(6)]):
    z_s = 650
    
    ## HRM
    asym3 = np.deg2rad(alphaAsym)
    crystal1 = optics.Crystal('C1', hkl=hkl, length=l_crystal[4], width=w_crystal[4],
                              z=z_s+310, E0=E0, alphaAsym=0, orientation=0, pol='s',
                              shapeError=shapeErrors[4])
    im_after_C1 = optics.PPM('im_after_C1', z=crystal1.z+.01, FOV=5e-3, N=256)

    crystal2 = optics.Crystal('C2', hkl=hkl, length=l_crystal[5], width=w_crystal[5],
                              z=crystal1.z+.05, E0=E0,alphaAsym=asym3, orientation=2, pol='s',
                              shapeError=shapeErrors[5])
    im_after_C2 = optics.PPM('im_after_C2', z=crystal2.z+.01, FOV=5e-3, N=256)
    
    im_before_MIR1 = optics.PPM('im_before_MIR1', z=crystal2.z + f1 - .1, FOV=2e-3, N=256)
    mir1 = optics.CurvedMirror('mir1', z=crystal2.z+f1, p=1e5, q=f2, length=1.0, width=5e-3, alpha=3.6e-3, orientation=0)
    im_after_MIR1 = optics.PPM('im_after_MIR1', z=mir1.z+.1, FOV=5e-3, N=256)

    # slit at focus
    slit = optics.Slit('Slit', z=mir1.z+f2, x_width=slit_width, y_width=2e-3)
    print('slit width: {} um'.format(slit.x_width*1e6))
    im_focus = optics.PPM('im_focus', z=mir1.z+f2 + 1e-3, FOV=50e-6, N=256)

    im_before_MIR2 = optics.PPM('im_before_MIR2', z=mir1.z+2*f2 -.1, FOV=2e-3, N=256)
    mir2 = optics.CurvedMirror('mir2', z=mir1.z+2*f2, p=f2, q=1e5, length=1.0, width=5e-3, alpha=3.6e-3, orientation=2)
    im_after_MIR2 = optics.PPM('im_after_MIR2', z=mir2.z+.1, FOV=2e-3, N=256)
    
    crystal3 = optics.Crystal('C3', hkl=hkl, length=1e-1, width=5e-3,
                              z=mir2.z+f1, E0=E0,alphaAsym=-asym3, orientation=2, pol='s')
    im_after_C3 = optics.PPM('im_after_C3', z=crystal3.z+.01, FOV=5e-3, N=256)

    crystal4 = optics.Crystal('C4', hkl=hkl, length=1e-1, width=5e-3,
                              z=crystal3.z+.05, E0=E0,alphaAsym=0, orientation=0, pol='s')
    im_out = optics.PPM('im_out', z=crystal4.z+.1, FOV=5e-3, N=256)
    
    HRM_devices = [crystal1,im_after_C1, crystal2,im_after_C2, im_before_MIR1,mir1,im_after_MIR1, slit,im_focus,
                   im_before_MIR2,mir2,im_after_MIR2, crystal3,im_after_C3, crystal4,im_out]
    
    return HRM_devices

def change_delta(devices, delta, crystal):
    for device in devices:
        if device.name == 'c{}'.format(crystal):
            device.delta = delta

def shift_z(mono_beamline, shift, oe):
    for device in mono_beamline.device_list:
        if device.name == oe:
            device.z += shift

def change_miscut(devices, eta_err, crystal):
    for i, device in enumerate(devices):
        if device.name == 'c{}'.format(crystal):
            devices[i] = optics.Crystal(device.name, hkl=device.hkl, length=device.length, width=device.width,
                                        z=device.z, E0=device.E0, alphaAsym=device.alphaAsym+eta_err,
                                        orientation=device.orientation, pol=device.pol, delta=device.delta)

def add_shapeError(devices, shapeError, crystal):
    for i, device in enumerate(devices):
        if device.name == 'c{}'.format(crystal):
            devices[i] = optics.Crystal(device.name, hkl=device.hkl, length=device.length, width=device.width,
                                        z=device.z, E0=device.E0, alphaAsym=device.alphaAsym,
                                        orientation=device.orientation, pol=device.pol, delta=device.delta,
                                        shapeError = shapeError)

def lens_energyError(devices, E):
    for i, device in enumerate(devices):
        if device.name[:3] == 'crl':
            devices[i] = optics.CRL('mir1', z=device.z, E0=E, f=device.f, diameter=device.diameter)

''' get info '''
def print_oe_pos(oe):
    print('{}, x:{}, y:{}, z:{}'.format(oe.name, oe.global_x, oe.global_y, oe.z))
    return oe.global_x, oe.global_y, oe.z

def get_width(pulse, image_name):
    # minima and maxima of the field of view (in microns) for imshow extent
    minx = np.round(np.min(pulse.x[image_name]) * 1e6)
    maxx = np.round(np.max(pulse.x[image_name]) * 1e6)
    miny = np.round(np.min(pulse.y[image_name]) * 1e6)
    maxy = np.round(np.max(pulse.y[image_name]) * 1e6)

    # calculate the profile
    profile = np.sum(np.abs(pulse.energy_stacks[image_name])**2, axis=2)
    x_lineout = np.sum(profile, axis=0)
    y_lineout = np.sum(profile, axis=1)

    cx, cy, wx, wy, fwx_guess, fwy_guess = pulse.beam_analysis(pulse.x[image_name],pulse.y[image_name],
                                                                       x_lineout, y_lineout)

    return cx, cy, wx, wy, fwx_guess, fwy_guess

def get_pulse(pulse, image_name, x_pos=0, y_pos=0, shift=None):
    minx = np.round(np.min(pulse.x[image_name]) * 1e6)
    maxx = np.round(np.max(pulse.x[image_name]) * 1e6)
    miny = np.round(np.min(pulse.y[image_name]) * 1e6)
    maxy = np.round(np.max(pulse.y[image_name]) * 1e6)

    # get number of pixels
    M = pulse.x[image_name].size
    N = pulse.y[image_name].size

    # calculate pixel sizes (microns)
    dx = (maxx - minx) / M
    dy = (maxy - miny) / N

    # calculate indices for the desired location
    x_index = int((x_pos - minx) / dx)
    y_index = int((y_pos - miny) / dy)

    # calculate temporal intensity
    y_data = np.abs(pulse.time_stacks[image_name][y_index, x_index, :]) ** 2

    shift = -pulse.t_axis[np.argmax(y_data)]

    # coarse shift for fitting
    if shift is not None:
        y_data = np.roll(y_data, int(shift/pulse.deltaT))

    # get gaussian stats
    centroid, sx = Util.gaussian_stats(pulse.t_axis, y_data)
    fwhm = int(sx * 2.355)

    # gaussian fit
    gauss_plot = Util.fit_gaussian(pulse.t_axis, centroid, sx)

    # shift again using fit result
    shift = -centroid
    if shift is not None:
        y_data = np.roll(y_data, int(shift/pulse.deltaT))
        gauss_plot = np.roll(gauss_plot, int(shift/pulse.deltaT))
        
    # [fs], normalized intensity [simulated], [Gaussian Fit]
    return pulse.t_axis, y_data/np.max(y_data), gauss_plot

def get_spectrum(pulse, image_name, x_pos=0, y_pos=0, integrated=False):
    minx = np.round(np.min(pulse.x[image_name]) * 1e6)
    maxx = np.round(np.max(pulse.x[image_name]) * 1e6)
    miny = np.round(np.min(pulse.y[image_name]) * 1e6)
    maxy = np.round(np.max(pulse.y[image_name]) * 1e6)

    # get number of pixels
    M = pulse.x[image_name].size
    N = pulse.y[image_name].size

    # calculate pixel sizes (microns)
    dx = (maxx - minx) / M
    dy = (maxy - miny) / N

    # calculate indices for the desired location
    x_index = int((x_pos - minx) / dx)
    y_index = int((y_pos - miny) / dy)

    # calculate spectral intensity
    if integrated:
        y_data = np.sum(np.abs(pulse.energy_stacks[image_name])**2, axis=(0,1))
    else:
        y_data = np.abs(pulse.energy_stacks[image_name][y_index,x_index,:])**2

    # get gaussian stats
    centroid, sx = Util.gaussian_stats(pulse.energy, y_data)
    fwhm = sx * 2.355

    # gaussian fit to plot
    gauss_plot = Util.fit_gaussian(pulse.energy, centroid, sx)

    # change label depending on bandwidth
    if fwhm >= 1:
        width_label = '%.1f eV FWHM' % fwhm
    elif fwhm > 1e-3:
        width_label = '%.1f meV FHWM' % (fwhm * 1e3)
    else:
        width_label = u'%.1f \u03BCeV FWHM' % (fwhm * 1e6)
    
    # [eV], normalized intensity [simulated], [Gaussian Fit]
    return pulse.energy - pulse.E0, y_data/np.max(y_data), gauss_plot