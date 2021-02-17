#!/usr/bin/env python
from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility

import os
try:
    __IPYTHON__
    import sys
    del sys.argv[1:]
except:
    pass


import srwl_bl
import srwlib
import srwlpy
import srwl_uti_smp
import uti_io
import math

from time import *
from copy import *
from array import *
from uti_plot import *

from util_Matt import Util
import numpy as np
import matplotlib
matplotlib.use('Agg')   # allows plot without X11 forwarding
import matplotlib.pyplot as plt


####### get info from SRW wavefront
def get_axis_sp(_wfr):
    # get axis in space (x, y); this leaves the pulse in time domain
    srwlpy.SetRepresElecField(_wfr, 't')
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
    axis_t = np.linspace(mesh_temp.eStart, mesh_temp.eFin, mesh_temp.ne)
    return axis_t

def get_intensity(_wfr, domain='t', polarization='total'):
    if domain == 'f':
        axis_x, axis_y = get_axis_sp(_wfr)            # get spatial axis; pulse is now in time domain
        ec = np.mean(get_axis_ev(_wfr))               # get energy center; pulse is now in frequency domain
    elif domain == 't':
        ec = np.mean(get_axis_ev(_wfr))               # get energy center; pulse is now in frequency domain
        axis_x, axis_y = get_axis_sp(_wfr)            # get spatial axis; pulse is now in time domain
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
    axis_x, axis_y = get_axis_sp(_wfr)            # get spatial axis; pulse is now in time domain
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
def plot_spatial_from_wf(_wfr):
    # plot wavefront projection (y vs x)
    img = get_intensity(_wfr, domain='t').sum(axis=-1)
    axis_x, axis_y = get_axis_sp(_wfr)
    plt.imshow(img,cmap='jet',
        extent = [axis_x.min()*1e6, axis_x.max()*1e6, axis_y.max()*1e6, axis_y.min()*1e6])
    plt.colorbar()
    plt.xlabel(r'x ($\mu$m)')
    plt.ylabel(r'y ($\mu$m)')

def plot_tprofile_from_wf(_wfr, if_short=1):
    # plot temporal profile (intensity vs time), if_short then only plot slices with > 1% intensity
    aw, axis_t, int0 = get_tprofile(_wfr)
    if if_short == 1:
        axis_t = axis_t[aw]
        int0 = int0[aw]
    plt.plot(axis_t*1e15, int0)
    plt.xlabel('time (fs)')
    plt.ylabel('intensity (a.u.)')

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
                  extent = [axis_t.min()*1e15, axis_t.max()*1e15, axis_sp.max()*1e6, axis_sp.min()*1e6])
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
    if if_short == 1:
        axis_ev = axis_ev[aw]
        int0 = int0[aw]
    plt.plot(axis_ev, int0)
    plt.xlabel('photon energy (eV)')
    plt.ylabel('intensity (a.u.)')

####### Fit
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

    return fwhm

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