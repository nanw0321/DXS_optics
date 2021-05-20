#!/usr/bin/env python
# coding: utf-8

# In[1]:


###### load modules
from Module_diagnostic_functions import *
from Module_propagation_functions import *


# In[2]:


import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--bl', type=str, required=False,
    help='beamline type: HHLM, HRM, HHLM_HRM')
parser.add_argument('--var_t', type=float, required=False,
    help='variable: pulse_duration [fs]')
parser.add_argument('--err_E', type=float, required=False,
    help='error: photon energy [meV]')
parser.add_argument('--err_f0', type=float, required=False,
    help='error: collimation lens focal distance [m]')
parser.add_argument('--err_f1', type=float, required=False,
    help='error: crystal-lens distance [m]')
parser.add_argument('--C', type=str, required=False,
    help='crystal name [string]')
parser.add_argument('--err_delta', type=float, required=False,
    help='error: crystal alignment [nrad]')
parser.add_argument('--err_miscut', type=float, required=False,
    help='error: crystal miscut [nrad]')
parser.add_argument('--if_log', type=int, required=False,
    help='condition: log scale plots')
parser.add_argument('--if_slice', type=int, required=False,
    help='condition: y=0 slice plots')
parser.add_argument('--if_close', type=int, required=False,
    help='condition: closed slit')
parser.add_argument('--job', type=int, required=False,
    help='job number')

args = parser.parse_args()

bl_name = args.bl
var_t = args.var_t
err_E = args.err_E
err_f0 = args.err_f0
err_f1 = args.err_f1

cname = args.C
err_delta = args.err_delta
err_miscut = args.err_miscut

if_log = args.if_log
if_slice = args.if_slice
if_close = args.if_close

job_num = args.job

err_name = ''
err_val_name = ''

if bl_name is None: bl_name = 'HRM'

if var_t is None: var_t = 400

if err_E is None: err_E = 0.
elif err_E != 0. : err_name += '_E'; err_val_name += '_{}meV'.format(round(err_E,3))

if err_f0 is None: err_f0 = 0.
elif err_f0 != 0.: err_name += '_f0'; err_val_name += '_{}m'.format(round(err_f0,3))

if err_f1 is None: err_f1 = 0.
elif err_f1 != 0.: err_name += '_f1'; err_val_name += '_{}m'.format(round(err_f1,3))

if cname is None: cname = 'HHLM1'

if err_delta is None: err_delta = 0.
elif err_delta != 0.: err_name += '_delta_{}'.format(cname); err_val_name += '_{}nrad'.format(round(err_delta,3))

if err_miscut is None: err_miscut = 0.
elif err_miscut != 0.: err_name += '_miscut_{}'.format(cname); err_val_name += '_{}nrad'.format(round(err_miscut,3))

if if_log is None: if_log = 1
if if_slice is None: if_slice = 0
if if_close is None: if_close = 0
    
if job_num is None: job_num = 0

print('arguments: {}'.format([var_t, err_E, err_f0, err_f1, err_delta, err_miscut]))
print('error type: {}, error value: {}'.format(err_name, err_val_name))
print('if_log {}, if_close {}, if_slice {}'.format(if_log, if_close, if_slice))


# In[3]:


##### sampling parameters
sigT = var_t * 1e-15/2.355
E = 9481.
d_slit = 7e-6

t_res =sigT*2.355/10              # time sampling resolution [s]
ev_window = 4/t_res *1e-15        # total energy window [eV]
ev_res = min(ev_window/800, 1e-3) # energy sampling resolution [eV]


range_x = 5e-3; range_y = 5e-3
nx = 1024; ny = 8; nz = 2*round(ev_window/ev_res/2)

x_res = range_x/nx
y_res = range_y/ny
t_window = t_res*nz          # total time window [s]

pulseRange = int(t_window/sigT)
factor = -1 # factor = 0.5

print('nx, ny, nz: {}'.format([nx, ny, nz]))
print('x resolution/range: {}/ {}um'.format(round(x_res*1e6,2), round(range_x*1e6,2)))
print('y resolution/range: {}/ {}um'.format(round(y_res*1e6,2), round(range_y*1e6,2)))
print('time resolution/range: {}/ {}fs'.format(round(t_res*1e15,2), round(t_window*1e15,2)))
print('energy resolution/range: {}/ {}meV'.format(round(ev_res*1e3,2), round(ev_window*1e3,2)))


# In[4]:


##### general lens calculations
name_CRL = ['CRL0', 'MIR1', 'MIR2']
fCRL_list = [290, 10, 10]
nCRL_list = [1, 1, 1]

##### general crystal calculations
# Si 220
Si220 = srwlib.SRWLOptCryst(_d_sp=1.9201374688016222, _psi0r=-1.0873035035585694e-05, _psi0i=1.8438837339536554e-07,
                         _psi_hr=-6.610167427198717e-06, _psi_hi=1.7824173540780476e-07,
                         _psi_hbr=-6.610167427198717e-06, _psi_hbi=1.7824173540780476e-07,
                         _tc=0.01, _ang_as=0)

# Si 440
Si440 = srwlib.SRWLOptCryst(_d_sp=0.9600687344008111, _psi0r=-1.0873035035585694e-05, _psi0i=1.8438837339536554e-07,
                         _psi_hr=-4.181686438547451e-06, _psi_hi=1.6100412693351052e-07,
                         _psi_hbr=-4.181686438547451e-06, _psi_hbi=1.6100412693351052e-07,
                         _tc=0.01, _ang_as=0)

thetaB220 = Si220.get_ang_inc(_e=E)
thetaB440 = Si440.get_ang_inc(_e=E)


# In[5]:


##### beamline crystal calculations
''' input parameters '''
name_HHLM = ['HHLM1', 'HHLM2', 'HHLM3', 'HHLM4']
thetaB_list_HHLM = np.array([thetaB220, thetaB440, thetaB440, thetaB220])
ang_asym_list_HHLM = np.deg2rad([17, 0, 0, -17])
ang_dif_pl_list_HHLM = np.array([np.pi/2, -np.pi/2, np.pi/2, -np.pi/2])

name_HRM = ['C1', 'C2', 'C3', 'C4']
thetaB_list_HRM = np.array([thetaB440, thetaB440, thetaB440, thetaB440])
ang_asym_list_HRM = np.deg2rad([0, -29.5, 0, 29.5])
ang_dif_pl_list_HRM = np.array([np.pi/2, -np.pi/2, -np.pi/2, np.pi/2])
  
''' var_param '''
nv_list_HHLM = []; tv_list_HHLM = []; ang_in_HHLM = []
for i in range(len(name_HHLM)):
    thetaB = thetaB_list_HHLM[i]
    ang_asym = ang_asym_list_HHLM[i]
    ang_dif_pl = ang_dif_pl_list_HHLM[i]
    nv, tv, ang_in = calc_crystal_orientation(thetaB, ang_asym, ang_dif_pl)
    nv_list_HHLM.append(nv)
    tv_list_HHLM.append(tv)
    ang_in_HHLM.append(ang_in)

print('HHLM var_param\n',
      'nv: \n{}\n'.format(np.round(nv_list_HHLM,3)),
      'tv: \n{}\n'.format(np.round(tv_list_HHLM,3)),
      'incident angle: \n{}\n'.format(np.round(ang_in_HHLM,3)))

nv_list_HRM = []; tv_list_HRM = []; ang_in_HRM = []
for i in range(len(name_HRM)):
    thetaB = thetaB_list_HRM[i]
    ang_asym = ang_asym_list_HRM[i]
    ang_dif_pl = ang_dif_pl_list_HRM[i]
    nv, tv, ang_in = calc_crystal_orientation(thetaB, ang_asym, ang_dif_pl)
    nv_list_HRM.append(nv)
    tv_list_HRM.append(tv)
    ang_in_HRM.append(ang_in)
    
print('HRM var_param\n',
      'nv: \n{}\n'.format(np.round(nv_list_HRM,3)),
      'tv: \n{}\n'.format(np.round(tv_list_HRM,3)),
      'incident angle: \n{}\n'.format(np.round(ang_in_HRM,3)))

''' scale factors '''
b_factor_HHLM = calc_b_factor(thetaB_list_HHLM,ang_asym_list_HHLM).max()
x_scaling_HHLM = calc_scale_factor(b_factor_HHLM)
t_stretching_HHLM = calc_t_stretching(thetaB_list_HHLM,ang_asym_list_HHLM, range_x=range_x).max()
z_scaling_HHLM = np.max([np.round(t_stretching_HHLM/t_window),1.])

print('HHLM scaling\n',
      'b factor: {}\n'.format(round(b_factor_HHLM,2)),
      't stretching: {}fs\n'.format(round(t_stretching_HHLM*1e15,2)),
      'x scale: {}, z scale: {}\n'.format(x_scaling_HHLM, z_scaling_HHLM))

b_factor_HRM = calc_b_factor(thetaB_list_HRM,ang_asym_list_HRM).max()
x_scaling_HRM = calc_scale_factor(b_factor_HRM)
t_stretching_HRM = calc_t_stretching(thetaB_list_HRM,ang_asym_list_HRM, range_x=range_x).max()
z_scaling_HRM = np.max([np.round(2*t_stretching_HRM/t_window),1.])

print('HRM scaling\n',
      'b factor: {}\n'.format(round(b_factor_HRM,2)),
      't stretching: {}fs\n'.format(round(t_stretching_HRM*1e15,2)),
      'x scale: {}, z scale: {}\n'.format(x_scaling_HRM, z_scaling_HRM))

''' HHLM positions '''
zlist = np.array([0, 139, 361, 500])    # crystal position in lab frame [mm]

# beam direction calculation
for i in range(len(zlist)):
    # loop through crystals and calculates beam direction
    if i == 0:
        v_in = [0,0,1]
        v_out = calc_direction_output(v_in, thetaB_list_HHLM[i], ang_dif_pl = ang_dif_pl_list_HHLM[i])
        beam_dir_HHLM = np.array([v_out])
    else:
        v_in = beam_dir_HHLM[i-1]
        v_out = calc_direction_output(v_in, thetaB_list_HHLM[i], ang_dif_pl = ang_dif_pl_list_HHLM[i])
        beam_dir_HHLM = np.append(beam_dir_HHLM, [v_out], axis=0)

# beam path calculation
driftz = zlist[1:] - zlist[:-1]
drift_list = []
beam_pos = np.zeros([1,3])

for i in range(len(driftz)):
    beam_dir = beam_dir_HHLM[i]
    drift = driftz[i]/beam_dir[-1]
    drift_list = np.append(drift_list, drift)
    beam_pos = np.append(beam_pos, [beam_pos[i]+drift*beam_dir], axis=0)
drift_list = np.append(drift_list, 5000.0-np.sum(drift_list))/1e3
# visualisation
plt.figure()
plt.plot(beam_pos[:,-1], -beam_pos[:,0],label='pulse between crystals')
plt.plot((np.asarray([min(beam_pos[:,-1])-100,max(beam_pos[:,-1])+100])),[0,0],'--',label='incident beam direction')
plt.legend(); plt.xlabel('z (mm)'); plt.ylabel('x (mm)'); plt.axis('equal')


# In[6]:


##### introduce errors
if err_E != 0.: print('energy error added in var_param automatically')
    
if err_f0!= 0.:
    fCRL_list[0] += err_f0
    print('collimation lens focal distance error added to var_param')

if err_f1!= 0.: print('crystal-lens distance error added in var_param automatically')
    
if err_delta!= 0:
    for i in range(len(name_HHLM)):
        if cname == name_HHLM[i]:
            thetaB = thetaB_list_HHLM[i]
            ang_asym = ang_asym_list_HHLM[i]
            ang_dif_pl = ang_dif_pl_list_HHLM[i]
            nv, tv, ang_in = calc_crystal_orientation(thetaB, ang_asym, ang_dif_pl, delta=err_delta/1e9)
            nv_list_HHLM[i] = nv
            tv_list_HHLM[i] = tv
            ang_in_HHLM[i] = ang_in
    for i in range(len(name_HRM)):
        if cname == name_HRM[i]:
            thetaB = thetaB_list_HRM[i]
            ang_asym = ang_asym_list_HRM[i]
            ang_dif_pl = ang_dif_pl_list_HRM[i]
            nv, tv, ang_in = calc_crystal_orientation(thetaB, ang_asym, ang_dif_pl, delta=err_delta/1e9)
            nv_list_HRM[i] = nv
            tv_list_HRM[i] = tv
            ang_in_HRM[i] = ang_in
    print('crystal alignment error added to var_param')
    
if err_miscut != 0:
    for i in range(len(name_HHLM)):
        if cname == name_HHLM[i]:
            ang_asym_list_HHLM[i] += err_miscut/1e9
            thetaB = thetaB_list_HHLM[i]
            ang_asym = ang_asym_list_HHLM[i]
            ang_dif_pl = ang_dif_pl_list_HHLM[i]
            nv, tv, ang_in = calc_crystal_orientation(thetaB, ang_asym, ang_dif_pl)
            nv_list_HHLM[i] = nv
            tv_list_HHLM[i] = tv
            ang_in_HHLM[i] = ang_in
    for i in range(len(name_HRM)):
        if cname == name_HRM[i]:
            ang_asym_list_HRM[i] += err_miscut/1e9
            thetaB = thetaB_list_HRM[i]
            ang_asym = ang_asym_list_HRM[i]
            ang_dif_pl = ang_dif_pl_list_HRM[i]
            nv, tv, ang_in = calc_crystal_orientation(thetaB, ang_asym, ang_dif_pl)
            nv_list_HRM[i] = nv
            tv_list_HRM[i] = tv
            ang_in_HRM[i] = ang_in
    print('crystal miscut error added to var_param\n')


# In[7]:


varParam = srwl_bl.srwl_uti_ext_options([
    ['name', 's', 'full_beamline', 'simulation name'],

#---Data Folder
    ['fdir', 's', '', 'folder (directory) name for reading-in input and saving output data files'],


    ['gbm_x', 'f', 0.0, 'average horizontal coordinates of waist [m]'],
    ['gbm_y', 'f', 0.0, 'average vertical coordinates of waist [m]'],
    ['gbm_z', 'f', 0.0, 'average longitudinal coordinate of waist [m]'],
    ['gbm_xp', 'f', 0.0, 'average horizontal angle at waist [rad]'],
    ['gbm_yp', 'f', 0.0, 'average verical angle at waist [rad]'],
    ['gbm_ave', 'f', E+err_E/1e3, 'average photon energy [eV]'],
    ['gbm_pen', 'f', 0.001, 'energy per pulse [J]'],
    ['gbm_rep', 'f', 1, 'rep. rate [Hz]'],
    ['gbm_pol', 'f', 2, 'polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left'],
    ['gbm_sx', 'f', 9.787229999999999e-06, 'rms beam size vs horizontal position [m] at waist (for intensity)'],
    ['gbm_sy', 'f', 9.787229999999999e-06, 'rms beam size vs vertical position [m] at waist (for intensity)'],
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
    ['op_r', 'f', 290.0, 'longitudinal position of the first optical element [m]'],
    # Former appParam:
    ['rs_type', 's', 'g', 'source type, (u) idealized undulator, (t), tabulated undulator, (m) multipole, (g) gaussian beam'],

#---Beamline optics:
    # CRL: crl
    ['op_CRL_foc_plane', 'f', 1, 'focalPlane'],
    ['op_CRL_delta', 'f', 3.791135e-06, 'refractiveIndex'],
    ['op_CRL_atten_len', 'f', 0.008387, 'attenuationLength'],
    ['op_CRL_shape', 'f', 1, 'shape'],
    ['op_CRL_apert_h', 'f', range_x, 'horizontalApertureSize'],
    ['op_CRL_apert_v', 'f', range_y, 'verticalApertureSize'],
    ['op_CRL_r_min', 'f', rCRL(fCRL_list[0], nCRL_list[0]), 'tipRadius'],
    ['op_CRL_wall_thick', 'f', 5e-05, 'tipWallThickness'],
    ['op_CRL_x', 'f', 0.0, 'horizontalOffset'],
    ['op_CRL_y', 'f', 0.0, 'verticalOffset'],
    ['op_CRL_n', 'i', nCRL_list[0], 'numberOfLenses'],

    # CRL_HHLM1/C1: drift
    ['op_CRL_HHLM1_L', 'f', 5.0, 'length'],
    ['op_CRL_C1_L', 'f', 10.0, 'length'],
    
    # HHLM1: crystal
    ['op_HHLM1_hfn', 's', '', 'heightProfileFile'],
    ['op_HHLM1_dim', 's', 'x', 'orientation'],
    ['op_HHLM1_d_sp', 'f', Si220.dSp, 'dSpacing'],
    ['op_HHLM1_psi0r', 'f', Si220.psi0r, 'psi0r'],
    ['op_HHLM1_psi0i', 'f', Si220.psi0i, 'psi0i'],
    ['op_HHLM1_psiHr', 'f', Si220.psiHr, 'psiHr'],
    ['op_HHLM1_psiHi', 'f', Si220.psiHi, 'psiHi'],
    ['op_HHLM1_psiHBr', 'f', Si220.psiHbr, 'psiHBr'],
    ['op_HHLM1_psiHBi', 'f', Si220.psiHbi, 'psiHBi'],
    ['op_HHLM1_tc', 'f', Si220.tc, 'crystalThickness'],
    ['op_HHLM1_uc', 'f', Si220.uc, 'useCase'],
    ['op_HHLM1_ang_as', 'f', ang_asym_list_HHLM[0], 'asymmetryAngle'],
    ['op_HHLM1_nvx', 'f', nv_list_HHLM[0][0], 'nvx'],
    ['op_HHLM1_nvy', 'f', nv_list_HHLM[0][1], 'nvy'],
    ['op_HHLM1_nvz', 'f', nv_list_HHLM[0][2], 'nvz'],
    ['op_HHLM1_tvx', 'f', tv_list_HHLM[0][0], 'tvx'],
    ['op_HHLM1_tvy', 'f', tv_list_HHLM[0][1], 'tvy'],
    ['op_HHLM1_ang', 'f', ang_in_HHLM[0], 'grazingAngle'],
    ['op_HHLM1_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_HHLM1_energy', 'f', E, 'energy'],
    ['op_HHLM1_diffractionAngle', 'f', ang_dif_pl_list_HHLM[0], 'diffractionAngle'],

    # HHLM1_HHLM2: drift
    ['op_HHLM1_HHLM2_L', 'f', drift_list[0], 'length'],

    # HHLM2: crystal
    ['op_HHLM2_hfn', 's', '', 'heightProfileFile'],
    ['op_HHLM2_dim', 's', 'x', 'orientation'],
    ['op_HHLM2_d_sp', 'f', Si440.dSp, 'dSpacing'],
    ['op_HHLM2_psi0r', 'f', Si440.psi0r, 'psi0r'],
    ['op_HHLM2_psi0i', 'f', Si440.psi0i, 'psi0i'],
    ['op_HHLM2_psiHr', 'f', Si440.psiHr, 'psiHr'],
    ['op_HHLM2_psiHi', 'f', Si440.psiHi, 'psiHi'],
    ['op_HHLM2_psiHBr', 'f', Si440.psiHbr, 'psiHBr'],
    ['op_HHLM2_psiHBi', 'f', Si440.psiHbi, 'psiHBi'],
    ['op_HHLM2_tc', 'f', Si440.tc, 'crystalThickness'],
    ['op_HHLM2_uc', 'f', Si440.uc, 'useCase'],
    ['op_HHLM2_ang_as', 'f', ang_asym_list_HHLM[1], 'asymmetryAngle'],
    ['op_HHLM2_nvx', 'f', nv_list_HHLM[1][0], 'nvx'],
    ['op_HHLM2_nvy', 'f', nv_list_HHLM[1][1], 'nvy'],
    ['op_HHLM2_nvz', 'f', nv_list_HHLM[1][2], 'nvz'],
    ['op_HHLM2_tvx', 'f', tv_list_HHLM[1][0], 'tvx'],
    ['op_HHLM2_tvy', 'f', tv_list_HHLM[1][1], 'tvy'],
    ['op_HHLM2_ang', 'f', ang_in_HHLM[1], 'grazingAngle'],
    ['op_HHLM2_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_HHLM2_energy', 'f', E, 'energy'],
    ['op_HHLM2_diffractionAngle', 'f', ang_dif_pl_list_HHLM[1], 'diffractionAngle'],

    # HHLM2_HHLM3: drift
    ['op_HHLM2_HHLM3_L', 'f', drift_list[1], 'length'],

    # HHLM3: crystal
    ['op_HHLM3_hfn', 's', '', 'heightProfileFile'],
    ['op_HHLM3_dim', 's', 'x', 'orientation'],
    ['op_HHLM3_d_sp', 'f', Si440.dSp, 'dSpacing'],
    ['op_HHLM3_psi0r', 'f', Si440.psi0r, 'psi0r'],
    ['op_HHLM3_psi0i', 'f', Si440.psi0i, 'psi0i'],
    ['op_HHLM3_psiHr', 'f', Si440.psiHr, 'psiHr'],
    ['op_HHLM3_psiHi', 'f', Si440.psiHi, 'psiHi'],
    ['op_HHLM3_psiHBr', 'f', Si440.psiHbr, 'psiHBr'],
    ['op_HHLM3_psiHBi', 'f', Si440.psiHbi, 'psiHBi'],
    ['op_HHLM3_tc', 'f', Si440.tc, 'crystalThickness'],
    ['op_HHLM3_uc', 'f', Si440.uc, 'useCase'],
    ['op_HHLM3_ang_as', 'f', ang_asym_list_HHLM[2], 'asymmetryAngle'],
    ['op_HHLM3_nvx', 'f', nv_list_HHLM[2][0], 'nvx'],
    ['op_HHLM3_nvy', 'f', nv_list_HHLM[2][1], 'nvy'],
    ['op_HHLM3_nvz', 'f', nv_list_HHLM[2][2], 'nvz'],
    ['op_HHLM3_tvx', 'f', tv_list_HHLM[2][0], 'tvx'],
    ['op_HHLM3_tvy', 'f', tv_list_HHLM[2][1], 'tvy'],
    ['op_HHLM3_ang', 'f', ang_in_HHLM[2], 'grazingAngle'],
    ['op_HHLM3_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_HHLM3_energy', 'f', E, 'energy'],
    ['op_HHLM3_diffractionAngle', 'f', ang_dif_pl_list_HHLM[2], 'diffractionAngle'],

    # HHLM3_HHLM4: drift
    ['op_HHLM3_HHLM4_L', 'f', drift_list[2], 'length'],

    # HHLM4: crystal
    ['op_HHLM4_hfn', 's', '', 'heightProfileFile'],
    ['op_HHLM4_dim', 's', 'x', 'orientation'],
    ['op_HHLM4_d_sp', 'f', Si220.dSp, 'dSpacing'],
    ['op_HHLM4_psi0r', 'f', Si220.psi0r, 'psi0r'],
    ['op_HHLM4_psi0i', 'f', Si220.psi0i, 'psi0i'],
    ['op_HHLM4_psiHr', 'f', Si220.psiHr, 'psiHr'],
    ['op_HHLM4_psiHi', 'f', Si220.psiHi, 'psiHi'],
    ['op_HHLM4_psiHBr', 'f', Si220.psiHbr, 'psiHBr'],
    ['op_HHLM4_psiHBi', 'f', Si220.psiHbi, 'psiHBi'],
    ['op_HHLM4_tc', 'f', Si220.tc, 'crystalThickness'],
    ['op_HHLM4_uc', 'f', Si220.uc, 'useCase'],
    ['op_HHLM4_ang_as', 'f', ang_asym_list_HHLM[3], 'asymmetryAngle'],
    ['op_HHLM4_nvx', 'f', nv_list_HHLM[3][0], 'nvx'],
    ['op_HHLM4_nvy', 'f', nv_list_HHLM[3][1], 'nvy'],
    ['op_HHLM4_nvz', 'f', nv_list_HHLM[3][2], 'nvz'],
    ['op_HHLM4_tvx', 'f', tv_list_HHLM[3][0], 'tvx'],
    ['op_HHLM4_tvy', 'f', tv_list_HHLM[3][1], 'tvy'],
    ['op_HHLM4_ang', 'f', ang_in_HHLM[3], 'grazingAngle'],
    ['op_HHLM4_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_HHLM4_energy', 'f', E, 'energy'],
    ['op_HHLM4_diffractionAngle', 'f', ang_dif_pl_list_HHLM[3], 'diffractionAngle'],
    
    # HHLM4_C1: drift
    ['op_HHLM4_C1_L', 'f', drift_list[3], 'length'],
    
    # C1: crystal
    ['op_C1_hfn', 's', '', 'heightProfileFile'],
    ['op_C1_dim', 's', 'x', 'orientation'],
    ['op_C1_d_sp', 'f', Si440.dSp, 'dSpacing'],
    ['op_C1_psi0r', 'f', Si440.psi0r, 'psi0r'],
    ['op_C1_psi0i', 'f', Si440.psi0i, 'psi0i'],
    ['op_C1_psiHr', 'f', Si440.psiHr, 'psiHr'],
    ['op_C1_psiHi', 'f', Si440.psiHi, 'psiHi'],
    ['op_C1_psiHBr', 'f', Si440.psiHbr, 'psiHBr'],
    ['op_C1_psiHBi', 'f', Si440.psiHbi, 'psiHBi'],
    ['op_C1_tc', 'f', Si440.tc, 'crystalThickness'],
    ['op_C1_uc', 'f', Si440.uc, 'useCase'],
    ['op_C1_ang_as', 'f', ang_asym_list_HRM[0], 'asymmetryAngle'],
    ['op_C1_nvx', 'f', nv_list_HRM[0][0], 'nvx'],
    ['op_C1_nvy', 'f', nv_list_HRM[0][1], 'nvy'],
    ['op_C1_nvz', 'f', nv_list_HRM[0][2], 'nvz'],
    ['op_C1_tvx', 'f', tv_list_HRM[0][0], 'tvx'],
    ['op_C1_tvy', 'f', tv_list_HRM[0][1], 'tvy'],
    ['op_C1_ang', 'f', ang_in_HRM[0], 'grazingAngle'],
    ['op_C1_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_C1_energy', 'f', E, 'energy'],
    ['op_C1_diffractionAngle', 'f', ang_dif_pl_list_HRM[0], 'diffractionAngle'],

    # C1_C2: drift
    ['op_C1_C2_L', 'f', 0.2, 'length'],

    # C2: crystal
    ['op_C2_hfn', 's', '', 'heightProfileFile'],
    ['op_C2_dim', 's', 'x', 'orientation'],
    ['op_C2_d_sp', 'f', Si440.dSp, 'dSpacing'],
    ['op_C2_psi0r', 'f', Si440.psi0r, 'psi0r'],
    ['op_C2_psi0i', 'f', Si440.psi0i, 'psi0i'],
    ['op_C2_psiHr', 'f', Si440.psiHr, 'psiHr'],
    ['op_C2_psiHi', 'f', Si440.psiHi, 'psiHi'],
    ['op_C2_psiHBr', 'f', Si440.psiHbr, 'psiHBr'],
    ['op_C2_psiHBi', 'f', Si440.psiHbi, 'psiHBi'],
    ['op_C2_tc', 'f', Si440.tc, 'crystalThickness'],
    ['op_C2_uc', 'f', Si440.uc, 'useCase'],
    ['op_C2_ang_as', 'f', ang_asym_list_HRM[1], 'asymmetryAngle'],
    ['op_C2_nvx', 'f', nv_list_HRM[1][0], 'nvx'],
    ['op_C2_nvy', 'f', nv_list_HRM[1][1], 'nvy'],
    ['op_C2_nvz', 'f', nv_list_HRM[1][2], 'nvz'],
    ['op_C2_tvx', 'f', tv_list_HRM[1][0], 'tvx'],
    ['op_C2_tvy', 'f', tv_list_HRM[1][1], 'tvy'],
    ['op_C2_ang', 'f', ang_in_HRM[1], 'grazingAngle'],
    ['op_C2_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_C2_energy', 'f', E, 'energy'],
    ['op_C2_diffractionAngle', 'f', ang_dif_pl_list_HRM[1], 'diffractionAngle'],

    # C2_MIR1: drift
    ['op_C2_MIR1_L', 'f', fCRL_list[1]+err_f1, 'length'],

    # MIR1: ellipsoidMirror
    ['op_MIR1_hfn', 's', 'None', 'heightProfileFile'],
    ['op_MIR1_dim', 's', 'x', 'orientation'],
    ['op_MIR1_p', 'f', 1000000000000.0, 'firstFocusLength'],
    ['op_MIR1_q', 'f', 10.0, 'focalLength'],
    ['op_MIR1_ang', 'f', 0.0036, 'grazingAngle'],
    ['op_MIR1_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_MIR1_size_tang', 'f', 0.01, 'tangentialSize'],
    ['op_MIR1_size_sag', 'f', 0.01, 'sagittalSize'],
    ['op_MIR1_nvx', 'f', 0.9999935200069984, 'normalVectorX'],
    ['op_MIR1_nvy', 'f', 0.0, 'normalVectorY'],
    ['op_MIR1_nvz', 'f', -0.0035999922240050387, 'normalVectorZ'],
    ['op_MIR1_tvx', 'f', 0.0035999922240050387, 'tangentialVectorX'],
    ['op_MIR1_tvy', 'f', 0.0, 'tangentialVectorY'],
    ['op_MIR1_x', 'f', 0.0, 'horizontalOffset'],
    ['op_MIR1_y', 'f', 0.0, 'verticalOffset'],

    # MIR1_Slit: drift
    ['op_MIR1_Slit_L', 'f', fCRL_list[1], 'length'],

    # Slit: aperture
    ['op_Slit_shape', 's', 'r', 'shape'],
    ['op_Slit_Dx', 'f', d_slit, 'horizontalSize'],
    ['op_Slit_Dy', 'f', 0.01, 'verticalSize'],
    ['op_Slit_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Slit_y', 'f', 0.0, 'verticalOffset'],

    # Slit_MIR2: drift
    ['op_Slit_MIR2_L', 'f', fCRL_list[2], 'length'],

    # MIR2: ellipsoidMirror
    ['op_MIR2_hfn', 's', '', 'heightProfileFile'],
    ['op_MIR2_dim', 's', 'x', 'orientation'],
    ['op_MIR2_p', 'f', 10.0, 'firstFocusLength'],
    ['op_MIR2_q', 'f', 1000000000000.0, 'focalLength'],
    ['op_MIR2_ang', 'f', 0.0036, 'grazingAngle'],
    ['op_MIR2_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_MIR2_size_tang', 'f', 0.01, 'tangentialSize'],
    ['op_MIR2_size_sag', 'f', 0.01, 'sagittalSize'],
    ['op_MIR2_nvx', 'f', 0.9999935200069984, 'normalVectorX'],
    ['op_MIR2_nvy', 'f', 0.0, 'normalVectorY'],
    ['op_MIR2_nvz', 'f', -0.0035999922240050387, 'normalVectorZ'],
    ['op_MIR2_tvx', 'f', -0.0035999922240050387, 'tangentialVectorX'],
    ['op_MIR2_tvy', 'f', 0.0, 'tangentialVectorY'],
    ['op_MIR2_x', 'f', 0.0, 'horizontalOffset'],
    ['op_MIR2_y', 'f', 0.0, 'verticalOffset'],

    # MIR2_C3: drift
    ['op_MIR2_C3_L', 'f', fCRL_list[2]+err_f1, 'length'],
    
    # C3: crystal
    ['op_C3_hfn', 's', '', 'heightProfileFile'],
    ['op_C3_dim', 's', 'x', 'orientation'],
    ['op_C3_d_sp', 'f', Si440.dSp, 'dSpacing'],
    ['op_C3_psi0r', 'f', Si440.psi0r, 'psi0r'],
    ['op_C3_psi0i', 'f', Si440.psi0i, 'psi0i'],
    ['op_C3_psiHr', 'f', Si440.psiHr, 'psiHr'],
    ['op_C3_psiHi', 'f', Si440.psiHi, 'psiHi'],
    ['op_C3_psiHBr', 'f', Si440.psiHbr, 'psiHBr'],
    ['op_C3_psiHBi', 'f', Si440.psiHbi, 'psiHBi'],
    ['op_C3_tc', 'f', Si440.tc, 'crystalThickness'],
    ['op_C3_uc', 'f', Si440.uc, 'useCase'],
    ['op_C3_ang_as', 'f', ang_asym_list_HRM[2], 'asymmetryAngle'],
    ['op_C3_nvx', 'f', nv_list_HRM[2][0], 'nvx'],
    ['op_C3_nvy', 'f', nv_list_HRM[2][1], 'nvy'],
    ['op_C3_nvz', 'f', nv_list_HRM[2][2], 'nvz'],
    ['op_C3_tvx', 'f', tv_list_HRM[2][0], 'tvx'],
    ['op_C3_tvy', 'f', tv_list_HRM[2][1], 'tvy'],
    ['op_C3_ang', 'f', ang_in_HRM[2], 'grazingAngle'],
    ['op_C3_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_C3_energy', 'f', E, 'energy'],
    ['op_C3_diffractionAngle', 'f', ang_dif_pl_list_HRM[2], 'diffractionAngle'],

    # C3_C4: drift
    ['op_C3_C4_L', 'f', 0.2, 'length'],

    # C4: crystal
    ['op_C4_hfn', 's', '', 'heightProfileFile'],
    ['op_C4_dim', 's', 'x', 'orientation'],
    ['op_C4_d_sp', 'f', Si440.dSp, 'dSpacing'],
    ['op_C4_psi0r', 'f', Si440.psi0r, 'psi0r'],
    ['op_C4_psi0i', 'f', Si440.psi0i, 'psi0i'],
    ['op_C4_psiHr', 'f', Si440.psiHr, 'psiHr'],
    ['op_C4_psiHi', 'f', Si440.psiHi, 'psiHi'],
    ['op_C4_psiHBr', 'f', Si440.psiHbr, 'psiHBr'],
    ['op_C4_psiHBi', 'f', Si440.psiHbi, 'psiHBi'],
    ['op_C4_tc', 'f', Si440.tc, 'crystalThickness'],
    ['op_C4_uc', 'f', Si440.uc, 'useCase'],
    ['op_C4_ang_as', 'f', ang_asym_list_HRM[3], 'asymmetryAngle'],
    ['op_C4_nvx', 'f', nv_list_HRM[3][0], 'nvx'],
    ['op_C4_nvy', 'f', nv_list_HRM[3][1], 'nvy'],
    ['op_C4_nvz', 'f', nv_list_HRM[3][2], 'nvz'],
    ['op_C4_tvx', 'f', tv_list_HRM[3][0], 'tvx'],
    ['op_C4_tvy', 'f', tv_list_HRM[3][1], 'tvy'],
    ['op_C4_ang', 'f', ang_in_HRM[3], 'grazingAngle'],
    ['op_C4_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_C4_energy', 'f', E, 'energy'],
    ['op_C4_diffractionAngle', 'f', ang_dif_pl_list_HRM[3], 'diffractionAngle'],

    # C4_After_HRM: drift
    ['op_C4_After_HRM_L', 'f', 10.0, 'length'],
    
#---Propagation parameters
#                                   [0][1] [2] [3][4] [5]  [6]  [7]  [8]  [9] [10] [11]
    ['op_CRL_pp', 'f',              [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'CRL'],
    ['op_CRL_HHLM1_pp', 'f',        [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'CRL_HHLM1'],
    ['op_CRL_C1_pp', 'f',           [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'CRL_C1'],
    ['op_HHLM1_pp', 'f',            [0, 0, 1.0, 0, 0, 1.0, x_scaling_HHLM, 1.0, 1.0, 0.0, 0.0, 0.0], 'HHLM1'],
    ['op_HHLM1_HHLM2_pp', 'f',      [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'HHLM1_HHLM2'],
    ['op_HHLM2_pp', 'f',            [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'HHLM2'],
    ['op_HHLM2_HHLM3_pp', 'f',      [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'HHLM2_HHLM3'],
    ['op_HHLM3_pp', 'f',            [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'HHLM3'],
    ['op_HHLM3_HHLM4_pp', 'f',      [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'HHLM3_HHLM4'],
    ['op_HHLM4_pp', 'f',            [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'HHLM4'],
    ['op_HHLM4_C1_pp', 'f',         [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'HHLM4_C1'],
    ['op_C1_pp', 'f',               [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'C1'],
    ['op_C1_C2_pp', 'f',            [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'C1_C2'],
    ['op_C2_pp', 'f',               [0, 0, 1.0, 0, 0, 8*x_scaling_HRM, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'C2'],
    ['op_C2_MIR1_pp', 'f',          [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'C2_MIR1'],
    ['op_MIR1_pp', 'f',             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'MIR1'],
    ['op_MIR1_Slit_pp', 'f',        [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'MIR1_Slit'],
    ['op_Slit_pp', 'f',             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'Slit'],
    ['op_Slit_MIR2_pp', 'f',        [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'Slit_MIR2'],
    ['op_MIR2_pp', 'f',             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'MIR2'],
    ['op_MIR2_C3_pp', 'f',          [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'MIR2_C3'],
    ['op_C3_pp', 'f',               [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'C3'],
    ['op_C3_C4_pp', 'f',            [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'C3_C4'],
    ['op_C4_pp', 'f',               [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'C4'],
    ['op_C4_After_HRM_pp', 'f',     [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], 'C4_After_HRM'],
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


# In[9]:


####### I/O
import pickle, os, shutil

## output directories
dir_output = 'output/'
dir_bl_out = dir_output+'{}/'.format(bl_name)    # HHLM, HRM, HHLM_HRM
dir_dur_out = dir_bl_out+'{}fs{}/'.format(round(sigT*2.355*1e15,1), err_name)    # pulse duration (100.0fs/)
dir_plot = dir_dur_out+'{}fs_{}meV_job{}{}/'.format(round(t_window*1e15,1), round(ev_window*1e3,1), job_num, err_val_name)    # sampling range + error value

for dir_name in [dir_output, dir_bl_out, dir_dur_out, dir_plot]:
    mkdir(dir_name); # print('directory created: {}'.format(dir_name))

print('job {}\n    directories created'.format(job_num))

## store varParam
varParam_name = 'varParam{}.pkl'.format(job_num)
with open(varParam_name, 'wb') as f:
    pickle.dump(varParam, f)

print('    varParam generation completed')

## script variables
if if_close == 0: i_start = 0
# i_start = 100

if bl_name == 'HHLM':
    z_scaling = z_scaling_HHLM
    main_str = "main_HHLM('varParam{}.pkl', z_scaling={}, dir_plot='{}', if_log={}, if_slice={},  i_start={})".format(
                                    job_num, z_scaling, dir_plot, if_log, if_slice, i_start)
if bl_name == 'HRM':
    z_scaling = z_scaling_HRM
    if if_close == 1: i_start = 5
    i_start = 100
    main_str = "main_HRM('varParam{}.pkl', z_scaling={}, if_close={}, dir_plot='{}', if_log={}, if_slice={}, i_start={})".format(
                                    job_num, z_scaling, if_close, dir_plot, if_log, if_slice, i_start)
if bl_name == 'HHLM_HRM':
    z_scaling = z_scaling_HHLM
    if if_close == 1: i_start = 9
    i_start = 100
    main_str = "main_HHLM_HRM('varParam{}.pkl', z_scaling={}, if_close={}, dir_plot='{}', if_log={}, if_slice={}, i_start={})".format(
                                    job_num, z_scaling, if_close, dir_plot, if_log, if_slice, i_start)

## gen script
script_name = 'script{}_{}.py'.format(job_num, if_close)

with open(script_name, 'w') as fh:
    fh.writelines("from Module_diagnostic_functions import *\n")
    fh.writelines("from Module_propagation_functions import *\n")
    fh.writelines("if __name__ == '__main__':\n")
    fh.writelines("    {}".format(main_str))

print('    python script generation completed')

## gen job
job_name = 'job{}_{}.slrm'.format(job_num, if_close)
dir_log = 'logs/'; mkdir(dir_log)

with open(job_name, 'w') as fh:
    fh.writelines("#!/bin/bash -l\n")
    fh.writelines("#SBATCH -N 1\n")
    fh.writelines("#SBATCH -t 00:30:00\n")
    fh.writelines("#SBATCH -q debug\n")
    fh.writelines("#SBATCH -L SCRATCH\n")
    fh.writelines("#SBATCH -C haswell\n")
    fh.writelines("#SBATCH -o {}{}fs_{}{}_{}.log\n".format(dir_log, var_t, bl_name, err_name+err_val_name,if_close))
    fh.writelines("#SBATCH --account=lcls\n")
    fh.writelines("export OMP_NUM_THREADS=64\n")
    fh.writelines("#export OMP_DISPLAY_ENV=true \n")
    fh.writelines("#export OMP_DISPLAY_AFFINITY=true \n")
    fh.writelines("module load python\n")
    fh.writelines("conda activate SRWmp\n")
    fh.writelines("srun -n 1 python {}\n".format(script_name))

print('    SLURM job script generation completed')

## submit job
os.system("sbatch {}".format(job_name)); print('job {}_{} submitted'.format(job_num, if_close))

