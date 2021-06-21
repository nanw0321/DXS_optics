import numpy as np
from array import array
from copy import deepcopy
from h5py import File
from genesis.writers import write_openpmd_wavefront_h5
from Module_diagnostic_functions import *

##### load wavefront
dir_wfr = '/global/cfs/cdirs/lcls/dxs_sase/'
DFL = dir_wfr+'stats_9.00_0.out.dfl'
name_wfr = 'wavefront.h5'

# Genesis parameters. These need to be known to populate the openPMD-wavefront metadata

HXREn=9 #keV
PARAM={'ncar': 251, 'dgrid': 0.0007, 'xlamds': 1240e-9/(HXREn*1000), 'zsep': 20, 'ntail':0, 'itdp':1}

from genesis.parsers import parse_genesis_dfl
D2 = parse_genesis_dfl(DFL, nx=PARAM['ncar'])

print(PARAM['xlamds']*PARAM['zsep'])

with File(name_wfr, 'w') as h5:
    write_openpmd_wavefront_h5(h5, dfl=D2, param=PARAM)

h5 = File(name_wfr, 'r')
dict(h5['data']['000000']['meshes']['electricField'].attrs)
print('SRW wavefront written')

#### To SRW
from pmd_srw import srw_wfr_from_openpmd_wavefront

with File(name_wfr, 'r') as h5:
    arEx, arEy, kwargs,  wfr_attrs = srw_wfr_from_openpmd_wavefront(h5['data']['000000']['meshes'],  iz_step=None)

wfr = srwlib.SRWLWfr(_arEx=arEx, _arEy=arEy, **kwargs)
wfr.__dict__.update(wfr_attrs)


##### mirror calculations
# variables
divergence_ratio = 2    # rough calculation using RMS beam width
mir_theta = 2.65e-3    # grazing angle [rad]
# mir_theta = np.deg2rad(30)

# calculations
beam_size = sigX*np.sqrt(1+(z0/zR)**2)    # RMS beamsize
divergence = 2*np.arctan(beam_size/2/z0)    # incident beam divergence

divergence_out = divergence_ratio * divergence    # output beam divergence
a_para = np.arctan((divergence_out - divergence)/4)*np.sin(mir_theta)/beam_size    # parabolic mirror coefficient, y = ax^2

z_mir1 = 815    # position along beamline
z_mir2 = (928-z_mir1)/np.cos(mir_theta)+z_mir1

d_short = beam_size/np.sin(mir_theta) / np.sin(divergence_out) * np.sin(mir_theta)
d_long = beam_size/np.sin(mir_theta) / np.sin(divergence_out) * np.sin(np.pi-divergence_out-mir_theta)

p_mir2 = (z_mir2-z_mir1) + (d_short + d_long)/2    # effective source distance for mirror 2
beam_size_out = p_mir2 * divergence_out

# mirror 2
a_para2 = np.arctan(divergence_out/3) / beam_size_out

mir_size_tan1 = beam_size*6/np.sin(mir_theta)    # tangential size of mirror (long end)
mir_size_tan2 = beam_size_out*6/np.sin(mir_theta)
mir_size_sag = 0.1    # sagittal size of mirror (short end)

mir_size1 = (np.floor(mir_size_tan1/1e-1)+1)*1e-1
mir_size2 = (np.floor(mir_size_tan2/1e-1)+1)*1e-1

# write mirror profile to file
x_mir1 = np.linspace(-mir_size1/2, mir_size1/2, np.int(mir_size1/2e-5))
x_mir2 = np.linspace(-mir_size2/2, mir_size2/2, np.int(mir_size2/2e-5))
h_mir1 = -a_para*np.square(x_mir1)
h_mir2 = a_para*np.square(x_mir2)


np.savetxt('convex_mir.dat', np.c_[x_mir1, h_mir1], delimiter="\t")
np.savetxt('concave_mir.dat', np.c_[x_mir2, h_mir2], delimiter="\t")


##### I/O
def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)

# I/O directories
dir_output = 'output/'; mkdir(dir_output)
dir_plot = dir_output+'Telescope/'; mkdir(dir_plot)


##### define optics
def set_optics_MIR1(v=None):
    el = []
    pp = []
    names = ['MIR1']
    for el_name in names:
        if el_name == 'MIR1':
            # MIR1: mirror 160.0m
            mirror_file = v.op_MIR1_hfn
            assert os.path.isfile(mirror_file), \
                'Missing input file {}, required by MIR1 beamline element'.format(mirror_file)
            el.append(srwlib.srwl_opt_setup_surf_height_1d(
                srwlib.srwl_uti_read_data_cols(mirror_file, "\t", 0, 1),
                _dim=v.op_MIR1_dim,
                _ang=abs(v.op_MIR1_ang),
                _amp_coef=v.op_MIR1_amp_coef,
                _size_x=v.op_MIR1_size_x,
                _size_y=v.op_MIR1_size_y,
            ))
            pp.append(v.op_MIR1_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_MIR1_MIR2(v=None):
    el = []
    pp = []
    names = ['MIR1_MIR2']
    for el_name in names:
        if el_name == 'MIR1_MIR2':
            # MIR1_MIR2: drift 160.0m
            el.append(srwlib.SRWLOptD(
                _L=v.op_MIR1_MIR2_L,
            ))
            pp.append(v.op_MIR1_MIR2_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_MIR2(v=None):
    el = []
    pp = []
    names = ['MIR2']
    for el_name in names:
        if el_name == 'MIR2':
            # MIR1: mirror 160.0m
            mirror_file = v.op_MIR2_hfn
            assert os.path.isfile(mirror_file), \
                'Missing input file {}, required by MIR1 beamline element'.format(mirror_file)
            el.append(srwlib.srwl_opt_setup_surf_height_1d(
                srwlib.srwl_uti_read_data_cols(mirror_file, "\t", 0, 1),
                _dim=v.op_MIR2_dim,
                _ang=abs(v.op_MIR2_ang),
                _amp_coef=v.op_MIR2_amp_coef,
                _size_x=v.op_MIR2_size_x,
                _size_y=v.op_MIR2_size_y,
            ))
            pp.append(v.op_MIR2_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_MIR2_Out(v=None):
    el = []
    pp = []
    names = ['MIR2_Out']
    for el_name in names:    
        if el_name == 'MIR2_Out':
            # MIR2_output: drift 273.0m
            el.append(srwlib.SRWLOptD(
                _L=v.op_MIR2_Out_L,
            ))
            pp.append(v.op_MIR2_Out_pp)
    pp.append(v.op_fin_pp)
    return srwlib.SRWLOptC(el, pp)


varParam = srwl_bl.srwl_uti_ext_options([
    ['name', 's', 'full_beamline', 'simulation name'],

#---Data Folder
    ['fdir', 's', '', 'folder (directory) name for reading-in input and saving output data files'],


    ['gbm_x', 'f', 0.0, 'average horizontal coordinates of waist [m]'],
    ['gbm_y', 'f', 0.0, 'average vertical coordinates of waist [m]'],
    ['gbm_z', 'f', 0.0, 'average longitudinal coordinate of waist [m]'],
    ['gbm_xp', 'f', 0.0, 'average horizontal angle at waist [rad]'],
    ['gbm_yp', 'f', 0.0, 'average verical angle at waist [rad]'],
    ['gbm_ave', 'f', E+err_E, 'average photon energy [eV]'],
    ['gbm_pen', 'f', 0.001, 'energy per pulse [J]'],
    ['gbm_rep', 'f', 1, 'rep. rate [Hz]'],
    ['gbm_pol', 'f', 2, 'polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left'],
    ['gbm_sx', 'f', sigX, 'rms beam size vs horizontal position [m] at waist (for intensity)'],
    ['gbm_sy', 'f', sigY, 'rms beam size vs vertical position [m] at waist (for intensity)'],
    ['gbm_st', 'f', sigT, 'rms pulse duration [s] (for intensity)'],
    ['gbm_mx', 'f', 0, 'transverse Gauss-Hermite mode order in horizontal direction'],
    ['gbm_my', 'f', 0, 'transverse Gauss-Hermite mode order in vertical direction'],
    ['gbm_ca', 's', 'c', 'treat _sigX, _sigY as sizes in [m] in coordinate representation (_presCA="c") or as angular divergences in [rad] in angular representation (_presCA="a")'],
    ['gbm_ft', 's', 't', 'treat _sigT as pulse duration in [s] in time domain/representation (_presFT="t") or as bandwidth in [eV] in frequency domain/representation (_presFT="f")'],

#---Calculation Types
    #Single-Electron Intensity distribution vs horizontal and vertical position
    ['si', '', '', 'calculate single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position', 'store_true'],
    #Single-Electron Wavefront Propagation
    ['ws', '', '', 'calculate single-electron (/ fully coherent) wavefront propagation', 'store_true'],
    #Multi-Electron (partially-coherent) Wavefront Propagation
    ['wm', '', '', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],

    ['w_e', 'f', -pulseRange*sigT/2, 'photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ef', 'f', pulseRange*sigT/2, 'final photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ne', 'i', nz, 'number of points vs photon energy for calculation of intensity distribution'],
    ['w_x', 'f', 0.0, 'central horizontal position [m] for calculation of intensity distribution'],
    ['w_rx', 'f', range_x, 'range of horizontal position [m] for calculation of intensity distribution'],
    ['w_nx', 'i', nx, 'number of points vs horizontal position for calculation of intensity distribution'],
    ['w_y', 'f', 0.0, 'central vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ry', 'f', range_y, 'range of vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ny', 'i', ny, 'number of points vs vertical position for calculation of intensity distribution'],
    ['w_smpf', 'f', factor, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_meth', 'i', 2, 'method to use for calculation of intensity distribution vs horizontal and vertical position: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['w_prec', 'f', 0.01, 'relative precision for calculation of intensity distribution vs horizontal and vertical position'],

    ['w_ft', 's', 't', 'presentation/domain: "f"- frequency (photon energy), "t"- time'],

    ['w_u', 'i', 2, 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['si_pol', 'i', 6, 'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['si_type', 'i', 0, 'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],
    
    ['si_fn', 's', 'res_int_se.dat', 'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],
    ['si_pl', 's', '', 'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
    ['ws_fni', 's', 'res_int_pr_se.dat', 'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],
    ['ws_pl', 's', '', 'plot the resulting intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    #to add options
    ['op_r', 'f', z0, 'longitudinal position of the first optical element [m]'],
    # Former appParam:
    ['rs_type', 's', 'g', 'source type, (u) idealized undulator, (t), tabulated undulator, (m) multipole, (g) gaussian beam'],

#---Beamline optics:
    # MIR1: mirror
    ['op_MIR1_hfn', 's', 'convex_mir.dat', 'heightProfileFile'],
    ['op_MIR1_dim', 's', 'x', 'orientation'],
    ['op_MIR1_ang', 'f', mir_theta, 'grazingAngle'],
    ['op_MIR1_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_MIR1_size_x', 'f', mir_size1, 'horizontalTransverseSize'],
    ['op_MIR1_size_y', 'f', range_y, 'verticalTransverseSize'],

    # MIR1_MIR2: drift
    ['op_MIR1_MIR2_L', 'f', z_mir2-z_mir1, 'length'],
    
    # MIR2: mirror
    ['op_MIR2_hfn', 's', 'concave_mir.dat', 'heightProfileFile'],
    ['op_MIR2_dim', 's', 'x', 'orientation'],
    ['op_MIR2_ang', 'f', mir_theta, 'grazingAngle'],
    ['op_MIR2_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_MIR2_size_x', 'f', mir_size2, 'horizontalTransverseSize'],
    ['op_MIR2_size_y', 'f', range_y, 'verticalTransverseSize'],

    # MIR2_Out: drift
    ['op_MIR2_Out_L', 'f', 8.0, 'length'],
    
#---Propagation parameters
#                                   [0][1] [2] [3][4] [5]  [6]  [7]  [8]  [9] [10] [11]
    ['op_MIR1_pp', 'f',             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'MIR1'],
    ['op_MIR1_MIR2_pp', 'f',        [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'MIR1_MIR2'],
    ['op_MIR2_pp', 'f',             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'MIR2'],
    ['op_MIR2_Out_pp', 'f',         [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'MIR2_Out'],
    ['op_fin_pp', 'f',              [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'final post-propagation (resize) parameters'],

    #[ 0]: Auto-Resize (1) or not (0) Before propagation
    #[ 1]: Auto-Resize (1) or not (0) After propagation
    #[ 2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
    #[ 3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
    #[ 4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
    #[ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
    #[ 6]: Horizontal Resolution modification factor at Resizing
    #[ 7]: Vertical Range modification factor at Resizing
    #[ 8]: Vertical Resolution modification factor at Resizing
    #[ 9]: Type of wavefront Shift before Resizing (not yet implemented)
    #[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
    #[11]: New Vertical wavefront Center position after Shift (not yet implemented)
    #[12]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Horizontal Coordinate
    #[13]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Vertical Coordinate
    #[14]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Longitudinal Coordinate
    #[15]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Horizontal Coordinate
    #[16]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Vertical Coordinate
])


#### Pseudo-voigt fit
def fit_pulse_position_pvoigt(_wfr):
    # Method to calculate the beam position
    axis_x, axis_y = get_axis_sp(_wfr)       # get spatial axis
    image = get_intensity(_wfr, domain='t').sum(axis=-1)
    projection_x = image.sum(axis=0)
    projection_y = image.sum(axis=1)
    
    mod = PseudoVoigtModel()
    pars_x = mod.guess(projection_x, x=axis_x)
    out_x = mod.fit(projection_x, pars_x, x=axis_x)
    out_x.plot()
    centroid_x = out_x.values['center']; fwhm_x = out_x.values['fwhm']; alpha_x = out_x.values['fraction']
    
    pars_y = mod.guess(projection_y, x=axis_y)
    out_y = mod.fit(projection_y, pars_y, x=axis_y)
    out_y.plot()
    centroid_y = out_y.values['center']; fwhm_y = out_y.values['fwhm']; alpha_y = out_y.values['fraction']
    return centroid_x, centroid_y, fwhm_x, fwhm_y, alpha_x, alpha_y

def fit_pulse_bandwidth_pvoigt(_wfr):
    # Method to calculate the bandwidth of the pulse
    aw, axis_ev, int0 = get_spectrum(_wfr)
    axis_ev = axis_ev[aw]; int0 = int0[aw]

    mod = PseudoVoigtModel()
    pars = mod.guess(int0, x=axis_ev)
    out = mod.fit(int0, pars, x=axis_ev)
    out.plot()
    centroid = out.values['center']; fwhm = out.values['fwhm']; alpha = out.values['fraction']

    return centroid, fwhm, alpha

def fit_pulse_duration_pvoigt(_wfr):
    # Method to calculate the temporal pulse structure
    aw, axis_t, int0 = get_tprofile(_wfr)    # the pulse is now in time domain
    axis_t = axis_t[aw] * 1e15; int0 = int0[aw]
    shift = int(int0.size/2 - int0.argmax()) # distance between array center and peak
    y_data = np.roll(int0, shift)

    mod = PseudoVoigtModel()
    pars = mod.guess(y_data, x=axis_t)
    out = mod.fit(y_data, pars, x=axis_t)
    out.plot()
    centroid = out.values['center']; fwhm = out.values['fwhm']; alpha = out.values['fraction']

    return centroid, fwhm, alpha


##### define diagnostic plot function
def plot_wfr_diagnostic(_wfr, label=None, dir_plot=None, if_slice=1, if_log=0, i=0, i_start=0):
    '''
    if_slice: y = 0 slice only or full projection
    if_log: log scale in imshow plots
    i: plot index, make I/O easier
    i_start: minimum plot index, only indices >= i_start will be plotted
    '''
    if i < i_start: print('plot skipped')
    else:
        np.seterr(divide = 'ignore')
        nx, ny, nz = get_dimension(_wfr)
        # I/O plot name
        nx, ny, nz = get_dimension(_wfr)    # get dimension
        if if_log == 1:
            pltname = '{}_{}_nx{}_ny{}_nz{}_log.png'.format(i,label,nx,ny,nz)
        else:
            pltname = '{}_{}_nx{}_ny{}_nz{}.png'.format(i,label,nx,ny,nz)

        ''' some calculations '''
        # space
        title_space = label + ', '
        axis_x, axis_y = get_axis_sp(_wfr)      # spatial axis [m]
        res_x = round(np.abs(axis_x[1] - axis_x[0])*1e6,2)  # spatial resolution
        try:
            cent_x, cent_y, fwhm_x, fwhm_y, alpha_x, alpha_y = fit_pulse_position_pvoigt(_wfr)   # fit for pulse position and width
            xstart = cent_x - fwhm_x*10; xfin = cent_x + fwhm_x*10
            ystart = cent_y - fwhm_y*10; yfin = cent_y + fwhm_y*10
            title_space+= 'x size: {}um\n{}% Lorentzian, '.format(round(fwhm_x*1e6,2), round(alpha_x*100,2))
        except:
            xstart = axis_x.min(); xfin = axis_x.max()
            ystart = axis_y.min(); yfin = axis_y.max()
        xstart = max(axis_x.min()*1e6, xstart*1e6); xfin = min(axis_x.max()*1e6, xfin*1e6)
        ystart = max(axis_y.min()*1e6, ystart*1e6); yfin = min(axis_y.max()*1e6, yfin*1e6)
        title_space += 'resolution: {}um'.format(res_x)

        # energy
        title_energy = ''
        axis_ev = get_axis_ev(_wfr)     # energy axis [eV]
        axis_ev = axis_ev - axis_ev[int(axis_ev.size/2)]        # shifted energy axis
        res_ev = round(np.abs(axis_ev[1] - axis_ev[0])*1e3,2)   # energy resolution
        try:
            cent_E, fwhm_E, alpha_E = fit_pulse_bandwidth_pvoigt(_wfr)      # fit for pulse central energy and bandwidth
            Estart = -fwhm_E*10*1e3; Efin = fwhm_E*10*1e3
            title_energy += 'bandwidth: {}meV\n{}% Lorentzian, '.format(round(fwhm_E*1e3,2), round(alpha_E*100,2))
        except:
            Estart = axis_ev.min()*1e3; Efin = axis_ev.max()*1e3
        Estart = max(axis_ev.min()*1e3, Estart); Efin = min(axis_ev.max()*1e3, Efin)
        title_energy += 'resolution: {}meV'.format(res_ev)

        # time
        title_time = ''
        axis_t = get_axis_t(_wfr)       # time axis [s]
        res_t = round(np.abs(axis_t[1] - axis_t[0])*1e15,2)     # time resolution
        try:
            cent_t, fwhm_t, alpha_t = fit_pulse_duration_pvoigt(_wfr)       # fit for pulse duration
            tstart = cent_t - fwhm_t*10; tfin = cent_t + fwhm_t*10
            title_time += 'duration: {}fs\n{}% Lorentzian, '.format(round(fwhm_t,2), round(alpha_t*100,2))
        except:
            tstart = axis_t.min()*1e15; tfin = axis_t.max()*1e15; fwhm_t = None
        tstart = max(axis_t.min()*1e15, tstart); tfin = min(axis_t.max()*1e15, tfin)
        try:
            tilt = round(fit_pulsefront_tilt(_wfr, dim='x'),2)     # [fs/um]
            title_time += 'tilt: {}fs/um, '.format(tilt)
        except:
            tilt = None
        title_time += 'resolution: {}fs'.format(res_t)

        # divergence
        try:
            h = 4.135667696e-15
            dif_lim = fwhm_t*1e-15 * fwhm_E/h*2
            print('del_t*del_E/h: {}'.format(round(dif_lim,5)))
        except:
            print('del_t*del_E/h cannot be calculated')
        
        t0 = time(); print('plotting {},\n'.format(label),
                           '[nx, ny, nz]: {}\n'.format([nx, ny, nz]),
                           'x range: {}um\n'.format(round((axis_x.max()-axis_x.min())*1e6,2)),
                           't range: {}fs\n'.format(round((axis_t.max()-axis_t.min())*1e15,2)),
                           'E range: {}meV\n'.format(round((axis_ev.max()-axis_ev.min())*1e3,2)))
        
        ''' plots '''
        plt.figure(figsize=(24,6))
        # space
        plt.subplot(1,3,1); plot_spatial_from_wf(_wfr, if_slice=if_slice, if_log=if_log)
#         plt.title(title_space); plt.axis('tight'); plt.xlim([-(range_x * x_scaling_HHLM/2)*1e6, (range_x * x_scaling_HHLM/2)*1e6])
#         plt.title(title_space); plt.axis('tight'); plt.xlim([-range_x/2*1e6, range_x/2*1e6])
        plt.title(title_space); plt.axis('tight'); plt.xlim([xstart, xfin])
        if if_slice != 1: plt.ylim([-range_y/2*1e6, range_y/2*1e6])

        # time
        plt.subplot(1,3,2); plot_tilt_from_wf(_wfr, ori='Horizontal', type='slice', if_log=if_log)
        plt.title(title_time); plt.xlim([tstart, tfin]); plt.ylim([xstart, xfin])

        # energy
        plt.subplot(1,3,3); plot_spatial_spectrum_from_wf(_wfr, ori='Horizontal', if_slice=0, if_log=if_log)
        plt.title(title_energy); plt.xlim([xstart, xfin]); plt.ylim([Estart, Efin])

        plt.savefig(dir_plot+pltname)
#         plt.close('all')
        np.seterr(divide = 'warn')
        print('plot lasted {}s\n'.format(round(time()-t0,2)))


##### main
def main(if_log=1, if_slice=1, i_start=0):
    # initialization
    v = srwl_bl.srwl_uti_parse_options(varParam, use_sys_argv=True)
    v.si = True
    srwl_bl.SRWLBeamline(_name=v.name).calc_all(v)

    #### incident beam
    wfr = srwlib.SRWLWfr(_arEx=arEx, _arEy=arEy, **kwargs)
    wfr.__dict__.update(wfr_attrs); i_plot = 1
    plot_wfr_diagnostic(wfr, label='input', dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    # Telescope
    label = 'after MIR1'; print('Propagating through {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_MIR1(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)
    
    label = 'before MIR2'; print('Propagating through {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_MIR1_MIR2(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'after MIR2'; print('Propagating through {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_MIR2(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)
    
    label = 'Output'; print('Propagating through {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_MIR2_Out(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    
##### execution
if __name__ == '__main__':
    time_stamp = time()
    main(if_log=0, if_slice=0, i_start=0)
    print('\n\neverything lasted: {}s'.format(round(time()-time_stamp,2)))
