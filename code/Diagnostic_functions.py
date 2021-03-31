#!/usr/bin/env python
from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility

import os
try:
    __IPYTHON__
    import sys
    del sys.argv[1:]
except:
    pass

from srwpy import srwl_bl
from srwpy import srwlib
from srwpy import srwlpy
from srwpy import srwl_uti_smp
from srwpy import uti_io
from srwpy.uti_plot import *
import math

from time import *
from copy import *
from array import *

from util_Matt import Util
import numpy as np
import matplotlib
# matplotlib.use('Agg')   # allows plot without X11 forwarding
import matplotlib.pyplot as plt

interpolation = 'antialiased'
# interpolation = 'none'

####### calculations
def calc_b_factor(thetaB, ang_asym):
    # calculates the b factor for crystal reflections [unit]
    b_factor = np.sin(thetaB+ang_asym)/ np.sin(thetaB-ang_asym)
    return b_factor

def calc_scale_factor(b_factor):
    # calculates the spatial scaling factor of crystal reflections [unit]
    scale_factors = np.array([1,2,3,4,5,6,8,9,10,12,15,16,18,20,24,25,27,30])
    if b_factor>1:
        idx = (np.abs(scale_factors - b_factor)).argmin()
    else:
        idx = (np.abs(scale_factors - 1/b_factor)).argmin()
    return scale_factors[idx]

def calc_t_stretching(thetaB, ang_asym, range_x=None):
    # calculates the stretched beam size in time after asymmetric reflection [s]
    t_stretching = np.abs(range_x / np.sin(thetaB-ang_asym) * (np.cos(thetaB-ang_asym)-np.cos(thetaB+ang_asym)) /3e8)
    return t_stretching

def calc_rot_mat_x(theta):
    # CW rotation if viewed along the x axis
    Rotation_matrix_x = np.asarray([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta), np.cos(theta)]
    ])
    return Rotation_matrix_x

def calc_rot_mat_z(theta):
    # CW rotation if viewed along the z axis
    Rotation_matrix_z = np.asarray([
        [np.cos(theta), -np.sin(theta),0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])
    return Rotation_matrix_z

def calc_rot_mat_btw_vecs(v1, v2):
    # calculates the rotation matrix that turns v1 to v2
    v1 = np.asarray(v1)
    v2 = np.asarray(v2)
    if np.any(v1==v2):
        Rotation_matrix = np.eye(3)
    else:
        a, b = (v1 / np.linalg.norm(v1)).reshape(3), (v2 / np.linalg.norm(v2)).reshape(3)
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        Rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return Rotation_matrix

def calc_crystal_orientation(thetaB, ang_asym, ang_dif_pl, delta=0.):
    '''calculates the crystal orientation, dif_pl defines the direction of reflected beam in the incident frame
        ang_dif_pl: [0 to +y, pi/2 to -x, pi to -y, -pi/2 to +x], in betweeen values are allowed
            Note: Oleg's "+x" axis is actually the "-x" axis, so pi/2 corresponds to horizontal reflection to the left
        delta: crystal alignment error, >0 turns the normal vector toward the reflected beam
        tv represents direction of crystal surface
        nv represents direction of crystal surface normal
    '''
    ang_incidence = thetaB-ang_asym
    nv = np.array([0, np.cos(ang_incidence), -np.sin(ang_incidence)])
    tv = np.array([0, np.sin(ang_incidence), np.cos(ang_incidence)])
    rot_mat = calc_rot_mat_x(delta)
    nv = np.dot(rot_mat, nv)
    tv = np.dot(rot_mat, tv)
    rot_mat = calc_rot_mat_z(ang_dif_pl)
    nv_rot = np.dot(rot_mat, nv)
    tv_rot = np.dot(rot_mat, tv)
    return nv_rot, tv_rot, ang_incidence

def calc_direction_output(v_in, thetaB, ang_dif_pl=0):
    # calculates the reflected beam direction in the lab frame
    nv, _, _ = calc_crystal_orientation(thetaB, 0, ang_dif_pl=ang_dif_pl)
    rot_mat = calc_rot_mat_btw_vecs(np.array([0,0,1]), v_in)
    nv_rot = np.dot(rot_mat, nv)    # crystal normal vector in the lab frame
    v_out = v_in - 2*np.dot(v_in, nv_rot)*nv_rot
    return v_out/np.linalg.norm(v_out)


####### get info from SRW wavefront
def get_dimension(_wfr):
    # get nx, ny, nz dimensions
    mesh_temp = deepcopy(_wfr.mesh)
    nx = mesh_temp.nx
    ny = mesh_temp.ny
    nz = mesh_temp.ne
    return nx, ny, nz

def get_axis_sp(_wfr):
    # get axis in space (x, y)
    mesh_temp = deepcopy(_wfr.mesh)
    axis_x = np.linspace(mesh_temp.xStart, mesh_temp.xFin, mesh_temp.nx)
    axis_y = np.linspace(mesh_temp.yStart, mesh_temp.yFin, mesh_temp.ny)
    return axis_x, axis_y

def get_axis_ev(_wfr):
    # get axis in energy (eV); this leaves the pulse in frequency domain
    srwlpy.SetRepresElecField(_wfr, 'f')
    mesh_temp = deepcopy(_wfr.mesh)
    axis_ev = np.linspace(mesh_temp.eStart, mesh_temp.eFin, mesh_temp.ne)
    return axis_ev

def get_axis_t(_wfr):
    # get axis in time (s); this leaves the pulse in time domain
    srwlpy.SetRepresElecField(_wfr, 't')
    mesh_temp = deepcopy(_wfr.mesh)
    axis_t = np.linspace(mesh_temp.eFin, mesh_temp.eStart, mesh_temp.ne)
    return axis_t

def get_intensity(_wfr, domain='t', polarization='total'):
    if domain == 'f':
        axis_x, axis_y = get_axis_sp(_wfr)            # get spatial axis
        ec = np.mean(get_axis_ev(_wfr))               # get energy center; pulse is now in frequency domain
        srwlpy.SetRepresElecField(_wfr, 'f')
    elif domain == 't':
        axis_x, axis_y = get_axis_sp(_wfr)            # get spatial axis
        ec = np.mean(get_axis_ev(_wfr))               # get energy center; pulse is now in frequency domain
        srwlpy.SetRepresElecField(_wfr, 't')
        
    xc = np.mean(axis_x); yc = np.mean(axis_y)    # get spatial centers

    # get 3d intensity [y:x:z] profile
    if polarization == 'total' or (polarization is None):
        pol = 6
    elif polarization == 'horizontal':
        pol = 0
    elif polarization == 'vertical':
        pol = 1
    else:
        raise ValueError(
            'unknown polarization value, shoule be "total" or "horizontal" or "vertical"')

    mesh_temp = deepcopy(_wfr.mesh)
    intensity = array('f', [0]*mesh_temp.ne*mesh_temp.nx*mesh_temp.ny)
    srwlpy.SetRepresElecField(_wfr, domain)           # just to be sure the pulse is in the right domain
    srwlpy.CalcIntFromElecField(intensity, _wfr, pol, 0, 6, ec, xc, yc)
    intensity = np.array(intensity, dtype='float32', copy=False)
    intensity.shape = (
        mesh_temp.ny, mesh_temp.nx, mesh_temp.ne)
    return intensity

def get_tprofile(_wfr):
    # get temporal profile (intensity vs time); this leaves the pulse in time domain
    axis_t = get_axis_t(_wfr)                     # get time axis; pulse is now in time domain
    int0 = get_intensity(_wfr, domain='t').sum(axis=(0,1))    # intensity per time slice

    aw = [a[0] for a in np.argwhere(int0>int0.max()*0.01)] # indicies for meaningful slices with intensity > 1% of the maximum
    aw = np.asarray(aw)
    return aw, axis_t, int0

def get_tilt(_wfr, ori='Vertical', type='sum'):
    # get the pulse front tilt in horizontal or vertical plane (y or x vs time)
    axis_x, axis_y = get_axis_sp(_wfr)            # get spatial axis
    intensity = get_intensity(_wfr, domain='t')
    if ori == 'Vertical':
        axis_sp = axis_y
        if type == 'sum':
            tilt = intensity.sum(axis=1)
        else:
            tilt = intensity[:,int(len(axis_x)/2),:]
    else:
        axis_sp = axis_x
        if type == 'sum':
            tilt = intensity.sum(axis=0)
        else:
            tilt = intensity[int(len(axis_y)/2),:,:]
    return axis_sp, tilt

def get_spectrum(_wfr):
    # get spectral info (intensity vs eV); this leaves the pulse in frequency domain
    axis_ev = get_axis_ev(_wfr)                   # get energy axis; pulse is now in frequency domain
    int0 = get_intensity(_wfr, domain='f').sum(axis=(0,1))    # intensity per energy slice

    aw = [a[0] for a in np.argwhere(int0>int0.max()*0.01)] # indicies for meaningful slices with intensity > 1% of the maximum
    aw = np.asarray(aw)
    return aw, axis_ev, int0


####### plot
def plot_spatial_from_wf(_wfr, if_slice=0, if_log=0):
    # plot wavefront projection (y vs x) or lineout (if_slice)
    nx, ny, nz = get_dimension(_wfr)
    img = get_intensity(_wfr, domain='t').sum(axis=-1)
    axis_x, axis_y = get_axis_sp(_wfr)
    if if_slice == 1:
        if nx >= ny:
            plt.plot(axis_x*1e6, img[int(ny/2)]/img[int(ny/2)].max())
            plt.xlabel(r'x ($\mu$m)')
        else:
            plt.plot(axis_y*1e6, img[:,int(nx/2)]/img[:,int(nx/2)].max())
            plt.xlabel(r'y ($\mu$m)')
        plt.ylabel('intensity (a.u.)')
    else:
        if if_log == 1:
            img = img/img.max()
            img = img + 1e-30
            img = np.log(img)
        plt.imshow(img,cmap='jet',
            extent = [axis_x.min()*1e6, axis_x.max()*1e6, axis_y.max()*1e6, axis_y.min()*1e6],
            interpolation=interpolation)
        plt.colorbar()
        if if_log == 1:
            cmin = np.max(img)-10
            plt.clim(cmin)
        plt.xlabel(r'x ($\mu$m)')
        plt.ylabel(r'y ($\mu$m)')

def plot_tprofile_from_wf(_wfr, if_short=1):
    # plot temporal profile (intensity vs time), if_short then only plot slices with > 1% intensity
    aw, axis_t, int0 = get_tprofile(_wfr)
    trange = (axis_t.max() - axis_t.min())*1e15
    npts = len(axis_t)

    if if_short == 1:
        axis_t = axis_t[aw]
        int0 = int0[aw]
    plt.plot(axis_t*1e15, int0/int0.max())
    plt.xlabel('time (fs)')
    plt.ylabel('intensity (a.u.)')
    plt.title('{} fs/{} pts'.format(round(trange,2), npts))

def plot_tilt_from_wf(_wfr, ori='Vertical', type='sum', if_log=0):
    # plot wavefront tilt (y or x vs time)
    axis_sp, tilt = get_tilt(_wfr, ori=ori, type=type)
    axis_t = get_axis_t(_wfr)

    def plot_tilt(axis_sp, tilt, axis_t, ori='Vertical', if_log=0):
        tilt = tilt/tilt.max()
        tilt = tilt + 1e-30
        if ori == 'Vertical': alabel = 'y'
        else: alabel = 'x'
        if if_log == 1:
            tilt = np.log(tilt)
        plt.imshow(tilt, cmap='jet',
                extent = [axis_t.max()*1e15, axis_t.min()*1e15, axis_sp.max()*1e6, axis_sp.min()*1e6],
                interpolation=interpolation)
        plt.colorbar()
        if if_log == 1:
            cmin = np.max(tilt)-10
            plt.clim(cmin)
        plt.axis('tight')
        plt.xlabel('time (fs)')
        plt.ylabel(alabel+r' ($\mu$m)')
    
    plot_tilt(axis_sp, tilt, axis_t, ori=ori, if_log=if_log)

def plot_spectrum_from_wf(_wfr, if_short=1):
    aw, axis_ev, int0 = get_spectrum(_wfr)
    ev_range = (axis_ev.max() - axis_ev.min())*1e3
    npts = len(axis_ev)
    ev_cent = axis_ev[int(len(axis_ev)/2)]
    
    if if_short == 1:
        axis_ev = axis_ev[aw]
        int0 = int0[aw]
    plt.plot( (axis_ev-ev_cent)*1e3, int0/int0.max())
    plt.xlabel('photon energy (meV) + {}eV'.format(ev_cent))
    plt.ylabel('intensity (a.u.)')
    plt.title('{} meV/{} pts'.format(round(ev_range,2), npts))

def plot_spatial_spectrum_from_wf(_wfr, ori='Vertical', if_slice=1, if_log=0):
    # spatial spectrum
    nx, ny, nz = get_dimension(_wfr)
    intensity_f = get_intensity(_wfr, domain='f')
    axis_ev = get_axis_ev(_wfr)
    ecent = axis_ev[int(axis_ev.size/2)]
    axis_x, axis_y = get_axis_sp(_wfr)
    
    if ori == 'Vertical':
        xlabel = 'y (um)'
        axis_sp = axis_y * 1e6
        if if_slice == 1:
            intensity_f = intensity_f[:,int(nx/2),:]
        else:
            intensity_f = intensity_f.sum(axis=1)
    elif ori == 'Horizontal':
        xlabel = 'x (um)'
        axis_sp = axis_x * 1e6
        if if_slice == 1:
            intensity_f = intensity_f[int(ny/2)]
        else:
            intensity_f = intensity_f.sum(axis=0)
    if if_slice == 1:
        title = 'spatial spectrum (central slice)'
    else:
        title = 'spatial spectrum (projection)'
    
    if if_log == 1:
        intensity_f = intensity_f/intensity_f.max()
        intensity_f = intensity_f + 1e-30
        intensity_f = np.log(intensity_f)
    plt.imshow(intensity_f, cmap='jet',
            extent=[axis_sp.min(), axis_sp.max(), (axis_ev.min()-ecent)*1e3, (axis_ev.max()-ecent)*1e3],
            interpolation=interpolation)
    plt.colorbar()
    if if_log == 1:
        cmin = np.max(intensity_f)-10
        plt.clim(cmin)
    plt.title(title); plt.xlabel(xlabel); plt.ylabel('meV around {}eV'.format(ecent))
    plt.axis('tight')

####### Fit
def fit_pulse_position(_wfr):
    # Method to calculate the beam position
    axis_x, axis_y = get_axis_sp(_wfr)       # get spatial axis
    image = get_intensity(_wfr, domain='t').sum(axis=-1)
    projection_x = image.sum(axis=0)
    projection_y = image.sum(axis=1)
    centroid_x, sigX = Util.gaussian_stats(axis_x, projection_x)
    centroid_y, sigY = Util.gaussian_stats(axis_y, projection_y)
    fwhm_x = sigX*2.355
    fwhm_y = sigY*2.355
    return centroid_x, centroid_y, fwhm_x, fwhm_y

def fit_pulse_duration(_wfr):
    # Method to calculate the temporal pulse structure
    aw, axis_t, int0 = get_tprofile(_wfr)    # the pulse is now in time domain
    axis_t = axis_t[aw] * 1e15; int0 = int0[aw]
    shift = int(int0.size/2 - int0.argmax()) # distance between array center and peak
    y_data = np.roll(int0, shift)

    # get Gaussian stats
    centroid, sigT = Util.gaussian_stats(axis_t, y_data)
    fwhm = int(sigT * 2.355)

    return centroid, fwhm

def fit_pulsefront_tilt(_wfr, dim='x'):
    # Method to calculate the pulse front tilt
    if dim == 'x':
        ori = 'Horizontal'
    else:
        ori = 'Vertical'
    axis_sp, tilt = get_tilt(_wfr, ori=ori, type='slice')
    axis_sp = axis_sp * 1e6
    axis_t = get_axis_t(_wfr) * 1e15

    # find peak at central spatial position
    index = np.argmax(tilt[int(axis_sp.size/2), :])

    # distance between array center and peak
    shift = int(np.size(tilt[int(axis_sp.size/2), :]) / 2 - index)

    profile = np.roll(tilt, shift, axis=1)

    # find peak (in time) at each position and put into fs units
    time_peaks = np.argmax(profile, axis=1) * (axis_t[1]-axis_t[0])

    # mask out anything outside the fwhm
    spatial_projection = np.sum(profile, axis=1)
    mask = spatial_projection>0.5*np.max(spatial_projection)

    time_peaks = time_peaks[mask]
    axis_sp = axis_sp[mask]

    # fit a line to the peaks
    p = np.polyfit(axis_sp, time_peaks, 1)
    # return slope (units are fs/micron)
    slope = p[0]

    return slope

def fit_pulse_bandwidth(_wfr):
    # Method to calculate the bandwidth of the pulse
    aw, axis_ev, int0 = get_spectrum(_wfr)
    axis_ev = axis_ev[aw]; int0 = int0[aw]

    # get gaussian stats
    centroid, sigE = Util.gaussian_stats(axis_ev, int0)
    fwhm = sigE * 2.355

    return centroid, fwhm

def fit_throughput(_wfr0, _wfr1):
    # Method to calculate the throughput at wfr1 relative to wfr0
    I0 = get_intensity(_wfr0, domain='t')
    I1 = get_intensity(_wfr1, domain='t')
    # this is a very rough treatment, for cross-comparison to be meaningful, need same ROI
    throughput = np.sum(I0)/ np.sum(I1)

    return throughput

def fit_central_energy(_wfr):
    # Method to calculate the central energy of a pulse
    aw, axis_ev, int0 = get_spectrum(_wfr)
    axis_ev = axis_ev[aw]; int0 = int0[aw]

    # get gaussian stats
    centroid, sigE = Util.gaussian_stats(axis_ev, int0)
    fwhm = sigE * 2.355

    return centroid