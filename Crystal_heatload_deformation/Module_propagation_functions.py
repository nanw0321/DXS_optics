#!/usr/bin/env python
# from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility

import os, pickle
try:
    __IPYTHON__
    import sys
    del sys.argv[1:]
except:
    pass

from Module_diagnostic_functions import *


####### I/O
def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)


####### define beamline optics
#### HHLM
def set_optics_CRL0_before_HHLM1(v=None):
    el = []
    pp = []
    names = ['CRL','CRL_HHLM1']
    for el_name in names:
        if el_name == 'CRL':
            # CRL: crl 290.0m
            el.append(srwlib.srwl_opt_setup_CRL(
                _foc_plane=v.op_CRL_foc_plane,
                _delta=v.op_CRL_delta,
                _atten_len=v.op_CRL_atten_len,
                _shape=v.op_CRL_shape,
                _apert_h=v.op_CRL_apert_h,
                _apert_v=v.op_CRL_apert_v,
                _r_min=v.op_CRL_r_min,
                _n=v.op_CRL_n,
                _wall_thick=v.op_CRL_wall_thick,
                _xc=v.op_CRL_x,
                _yc=v.op_CRL_y,
            ))
            pp.append(v.op_CRL_pp)
        elif el_name == 'CRL_HHLM1':
            # CRL_HHLM1: drift 290.0m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CRL_HHLM1_L,
            ))
            pp.append(v.op_CRL_HHLM1_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_HHLM1_before_HHLM2(v=None):
    el = []
    pp = []
    names = ['HHLM1', 'HHLM1_HHLM2']
    for el_name in names:
        if el_name == 'HHLM1':
            # HHLM1: crystal 295.0m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_HHLM1_d_sp,
                _psi0r=v.op_HHLM1_psi0r,
                _psi0i=v.op_HHLM1_psi0i,
                _psi_hr=v.op_HHLM1_psiHr,
                _psi_hi=v.op_HHLM1_psiHi,
                _psi_hbr=v.op_HHLM1_psiHBr,
                _psi_hbi=v.op_HHLM1_psiHBi,
                _tc=v.op_HHLM1_tc,
                _ang_as=v.op_HHLM1_ang_as,
                _nvx=v.op_HHLM1_nvx,
                _nvy=v.op_HHLM1_nvy,
                _nvz=v.op_HHLM1_nvz,
                _tvx=v.op_HHLM1_tvx,
                _tvy=v.op_HHLM1_tvy,
                _uc=v.op_HHLM1_uc,
                _e_avg=v.op_HHLM1_energy,
                _ang_roll=v.op_HHLM1_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_HHLM1_pp)
        elif el_name == 'HHLM1_HHLM2':
            # HHLM1_HHLM2: drift
            el.append(srwlib.SRWLOptD(
                _L=v.op_HHLM1_HHLM2_L,
            ))
            pp.append(v.op_HHLM1_HHLM2_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_HHLM2_before_HHLM3(v=None):
    el = []
    pp = []
    names = ['HHLM2', 'HHLM2_HHLM3']
    for el_name in names:
        if el_name == 'HHLM2':
            # HHLM2: crystal
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_HHLM2_d_sp,
                _psi0r=v.op_HHLM2_psi0r,
                _psi0i=v.op_HHLM2_psi0i,
                _psi_hr=v.op_HHLM2_psiHr,
                _psi_hi=v.op_HHLM2_psiHi,
                _psi_hbr=v.op_HHLM2_psiHBr,
                _psi_hbi=v.op_HHLM2_psiHBi,
                _tc=v.op_HHLM2_tc,
                _ang_as=v.op_HHLM2_ang_as,
                _nvx=v.op_HHLM2_nvx,
                _nvy=v.op_HHLM2_nvy,
                _nvz=v.op_HHLM2_nvz,
                _tvx=v.op_HHLM2_tvx,
                _tvy=v.op_HHLM2_tvy,
                _uc=v.op_HHLM2_uc,
                _e_avg=v.op_HHLM2_energy,
                _ang_roll=v.op_HHLM2_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_HHLM2_pp)
        elif el_name == 'HHLM2_HHLM3':
            # HHLM2_HHLM3: drift
            el.append(srwlib.SRWLOptD(
                _L=v.op_HHLM2_HHLM3_L,
            ))
            pp.append(v.op_HHLM2_HHLM3_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_HHLM3_before_HHLM4(v=None):
    el = []
    pp = []
    names = ['HHLM3', 'HHLM3_HHLM4']
    for el_name in names:
        if el_name == 'HHLM3':
            # HHLM3: crystal
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_HHLM3_d_sp,
                _psi0r=v.op_HHLM3_psi0r,
                _psi0i=v.op_HHLM3_psi0i,
                _psi_hr=v.op_HHLM3_psiHr,
                _psi_hi=v.op_HHLM3_psiHi,
                _psi_hbr=v.op_HHLM3_psiHBr,
                _psi_hbi=v.op_HHLM3_psiHBi,
                _tc=v.op_HHLM3_tc,
                _ang_as=v.op_HHLM3_ang_as,
                _nvx=v.op_HHLM3_nvx,
                _nvy=v.op_HHLM3_nvy,
                _nvz=v.op_HHLM3_nvz,
                _tvx=v.op_HHLM3_tvx,
                _tvy=v.op_HHLM3_tvy,
                _uc=v.op_HHLM3_uc,
                _e_avg=v.op_HHLM3_energy,
                _ang_roll=v.op_HHLM3_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_HHLM3_pp)
        elif el_name == 'HHLM3_HHLM4':
            # HHLM3_HHLM4: drift
            el.append(srwlib.SRWLOptD(
                _L=v.op_HHLM3_HHLM4_L,
            ))
            pp.append(v.op_HHLM3_HHLM4_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_HHLM4_output(v=None):
    el = []
    pp = []
    names = ['HHLM4']
    for el_name in names:
        if el_name == 'HHLM4':
            # HHLM4: crystal
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_HHLM4_d_sp,
                _psi0r=v.op_HHLM4_psi0r,
                _psi0i=v.op_HHLM4_psi0i,
                _psi_hr=v.op_HHLM4_psiHr,
                _psi_hi=v.op_HHLM4_psiHi,
                _psi_hbr=v.op_HHLM4_psiHBr,
                _psi_hbi=v.op_HHLM4_psiHBi,
                _tc=v.op_HHLM4_tc,
                _ang_as=v.op_HHLM4_ang_as,
                _nvx=v.op_HHLM4_nvx,
                _nvy=v.op_HHLM4_nvy,
                _nvz=v.op_HHLM4_nvz,
                _tvx=v.op_HHLM4_tvx,
                _tvy=v.op_HHLM4_tvy,
                _uc=v.op_HHLM4_uc,
                _e_avg=v.op_HHLM4_energy,
                _ang_roll=v.op_HHLM4_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_HHLM4_pp)
    pp.append(v.op_fin_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_HHLM4_before_C1(v=None):
    el = []
    pp = []
    names = ['HHLM4', 'HHLM4_C1']
    for el_name in names:
        if el_name == 'HHLM4':
            # HHLM4: crystal
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_HHLM4_d_sp,
                _psi0r=v.op_HHLM4_psi0r,
                _psi0i=v.op_HHLM4_psi0i,
                _psi_hr=v.op_HHLM4_psiHr,
                _psi_hi=v.op_HHLM4_psiHi,
                _psi_hbr=v.op_HHLM4_psiHBr,
                _psi_hbi=v.op_HHLM4_psiHBi,
                _tc=v.op_HHLM4_tc,
                _ang_as=v.op_HHLM4_ang_as,
                _nvx=v.op_HHLM4_nvx,
                _nvy=v.op_HHLM4_nvy,
                _nvz=v.op_HHLM4_nvz,
                _tvx=v.op_HHLM4_tvx,
                _tvy=v.op_HHLM4_tvy,
                _uc=v.op_HHLM4_uc,
                _e_avg=v.op_HHLM4_energy,
                _ang_roll=v.op_HHLM4_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_HHLM4_pp)
        elif el_name == 'HHLM4_C1':
            # HHLM3_HHLM4: drift
            el.append(srwlib.SRWLOptD(
                _L=v.op_HHLM4_C1_L,
            ))
            pp.append(v.op_HHLM4_C1_pp)
    return srwlib.SRWLOptC(el, pp)


#### HRM
## CC1
def set_optics_CRL0_before_C1(v=None):
    el = []
    pp = []
    names = ['CRL','CRL_C1']
    for el_name in names:
        if el_name == 'CRL':
            # CRL: crl 290.0m
            el.append(srwlib.srwl_opt_setup_CRL(
                _foc_plane=v.op_CRL_foc_plane,
                _delta=v.op_CRL_delta,
                _atten_len=v.op_CRL_atten_len,
                _shape=v.op_CRL_shape,
                _apert_h=v.op_CRL_apert_h,
                _apert_v=v.op_CRL_apert_v,
                _r_min=v.op_CRL_r_min,
                _n=v.op_CRL_n,
                _wall_thick=v.op_CRL_wall_thick,
                _xc=v.op_CRL_x,
                _yc=v.op_CRL_y,
            ))
            pp.append(v.op_CRL_pp)
        elif el_name == 'CRL_C1':
            # CRL_C1: drift 290.0m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CRL_C1_L,
            ))
            pp.append(v.op_CRL_C1_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_C1_before_C2(v=None):
    el = []
    pp = []
    names = ['C1', 'C1_C2']
    for el_name in names:
        if el_name == 'C1':
            # C1: crystal 300.0m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_C1_d_sp,
                _psi0r=v.op_C1_psi0r,
                _psi0i=v.op_C1_psi0i,
                _psi_hr=v.op_C1_psiHr,
                _psi_hi=v.op_C1_psiHi,
                _psi_hbr=v.op_C1_psiHBr,
                _psi_hbi=v.op_C1_psiHBi,
                _tc=v.op_C1_tc,
                _ang_as=v.op_C1_ang_as,
                _nvx=v.op_C1_nvx,
                _nvy=v.op_C1_nvy,
                _nvz=v.op_C1_nvz,
                _tvx=v.op_C1_tvx,
                _tvy=v.op_C1_tvy,
                _uc=v.op_C1_uc,
                _e_avg=v.op_C1_energy,
                _ang_roll=v.op_C1_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_C1_pp)

        elif el_name == 'C1_C2':
            # C1_C2: drift 300.0m
            el.append(srwlib.SRWLOptD(
                _L=v.op_C1_C2_L,
            ))
            pp.append(v.op_C1_C2_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_C2_before_CRL1(v=None):
    el = []
    pp = []
    names = ['C2', 'C2_CRL1']
    for el_name in names:
        if el_name == 'C2':
            # C2: crystal 300.2m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_C2_d_sp,
                _psi0r=v.op_C2_psi0r,
                _psi0i=v.op_C2_psi0i,
                _psi_hr=v.op_C2_psiHr,
                _psi_hi=v.op_C2_psiHi,
                _psi_hbr=v.op_C2_psiHBr,
                _psi_hbi=v.op_C2_psiHBi,
                _tc=v.op_C2_tc,
                _ang_as=v.op_C2_ang_as,
                _nvx=v.op_C2_nvx,
                _nvy=v.op_C2_nvy,
                _nvz=v.op_C2_nvz,
                _tvx=v.op_C2_tvx,
                _tvy=v.op_C2_tvy,
                _uc=v.op_C2_uc,
                _e_avg=v.op_C2_energy,
                _ang_roll=v.op_C2_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_C2_pp)
        elif el_name == 'C2_CRL1':
            # C2_CRL1: drift 300.2m
            el.append(srwlib.SRWLOptD(
                _L=v.op_C2_CRL1_L,
            ))
            pp.append(v.op_C2_CRL1_pp)
    return srwlib.SRWLOptC(el, pp)


## telescope
def set_optics_CRL1_focus(v=None):
    el = []
    pp = []
    names = ['CRL1', 'CRL1_Slit']
    for el_name in names:
        if el_name == 'CRL1':
            # CRL1: crl 310.2m
            el.append(srwlib.srwl_opt_setup_CRL(
                _foc_plane=v.op_CRL1_foc_plane,
                _delta=v.op_CRL1_delta,
                _atten_len=v.op_CRL1_atten_len,
                _shape=v.op_CRL1_shape,
                _apert_h=v.op_CRL1_apert_h,
                _apert_v=v.op_CRL1_apert_v,
                _r_min=v.op_CRL1_r_min,
                _n=v.op_CRL1_n,
                _wall_thick=v.op_CRL1_wall_thick,
                _xc=v.op_CRL1_x,
                _yc=v.op_CRL1_y,
            ))
            pp.append(v.op_CRL1_pp)
        elif el_name == 'CRL1_Slit':
            # CRL1_Slit: drift 310.2m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CRL1_Slit_L,
            ))
            pp.append(v.op_CRL1_Slit_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_slit(v=None):
    el = []
    pp = []
    names = ['Slit']
    for el_name in names:
        if el_name == 'Slit':
            # Slit: aperture 320.2m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_Slit_shape,
                _ap_or_ob='a',
                _Dx=v.op_Slit_Dx,
                _Dy=v.op_Slit_Dy,
                _x=v.op_Slit_x,
                _y=v.op_Slit_y,
            ))
            pp.append(v.op_Slit_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_slit_before_CRL2(v=None):
    el = []
    pp = []
    names = ['Slit_CRL2']
    for el_name in names:
        if el_name == 'Slit_CRL2':
            # Slit_CRL2: drift 320.2m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Slit_CRL2_L,
            ))
            pp.append(v.op_Slit_CRL2_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_CRL2_before_C3(v=None):
    el = []
    pp = []
    names = ['CRL2', 'CRL2_C3']
    for el_name in names:
        if el_name == 'CRL2':
            # CRL2: crl 330.2m
            el.append(srwlib.srwl_opt_setup_CRL(
                _foc_plane=v.op_CRL2_foc_plane,
                _delta=v.op_CRL2_delta,
                _atten_len=v.op_CRL2_atten_len,
                _shape=v.op_CRL2_shape,
                _apert_h=v.op_CRL2_apert_h,
                _apert_v=v.op_CRL2_apert_v,
                _r_min=v.op_CRL2_r_min,
                _n=v.op_CRL2_n,
                _wall_thick=v.op_CRL2_wall_thick,
                _xc=v.op_CRL2_x,
                _yc=v.op_CRL2_y,
            ))
            pp.append(v.op_CRL2_pp)
        elif el_name == 'CRL2_C3':
            # CRL2_C3: drift 330.2m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CRL2_C3_L,
            ))
            pp.append(v.op_CRL2_C3_pp)
    return srwlib.SRWLOptC(el, pp)


## CC2
def set_optics_C3_before_C4(v=None):
    el = []
    pp = []
    names = ['C3', 'C3_C4']
    for el_name in names:
        if el_name == 'C3':
            # C3: crystal 340.2m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_C3_d_sp,
                _psi0r=v.op_C3_psi0r,
                _psi0i=v.op_C3_psi0i,
                _psi_hr=v.op_C3_psiHr,
                _psi_hi=v.op_C3_psiHi,
                _psi_hbr=v.op_C3_psiHBr,
                _psi_hbi=v.op_C3_psiHBi,
                _tc=v.op_C3_tc,
                _ang_as=v.op_C3_ang_as,
                _nvx=v.op_C3_nvx,
                _nvy=v.op_C3_nvy,
                _nvz=v.op_C3_nvz,
                _tvx=v.op_C3_tvx,
                _tvy=v.op_C3_tvy,
                _uc=v.op_C3_uc,
                _e_avg=v.op_C3_energy,
                _ang_roll=v.op_C3_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_C3_pp)

        elif el_name == 'C3_C4':
            # C3_C4: drift 340.2m
            el.append(srwlib.SRWLOptD(
                _L=v.op_C3_C4_L,
            ))
            pp.append(v.op_C3_C4_pp)
    return srwlib.SRWLOptC(el, pp)


def set_optics_C4_output(v=None):
    el = []
    pp = []
    names = ['C4']
    for el_name in names:
        if el_name == 'C4':
            # C4: crystal 340.4m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_C4_d_sp,
                _psi0r=v.op_C4_psi0r,
                _psi0i=v.op_C4_psi0i,
                _psi_hr=v.op_C4_psiHr,
                _psi_hi=v.op_C4_psiHi,
                _psi_hbr=v.op_C4_psiHBr,
                _psi_hbi=v.op_C4_psiHBi,
                _tc=v.op_C4_tc,
                _ang_as=v.op_C4_ang_as,
                _nvx=v.op_C4_nvx,
                _nvy=v.op_C4_nvy,
                _nvz=v.op_C4_nvz,
                _tvx=v.op_C4_tvx,
                _tvy=v.op_C4_tvy,
                _uc=v.op_C4_uc,
                _e_avg=v.op_C4_energy,
                _ang_roll=v.op_C4_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_C4_pp)
    pp.append(v.op_fin_pp)
    return srwlib.SRWLOptC(el, pp)


####### propagation
#### HHLM
def main_HHLM(varParam_name, z_scaling=10.0, dir_plot=None, if_log=1, if_slice=1,  i_start=0):
    time_stamp=time()

    # initialization
    with open(varParam_name, 'rb') as f:
        varParam = pickle.load(f)
    v = srwl_bl.srwl_uti_parse_options(varParam, use_sys_argv=True)
    v.si = True
    srwl_bl.SRWLBeamline(_name=v.name).calc_all(v)

    # incident beam
    wfr = v.w_res; i_plot = 1
    plot_wfr_diagnostic(wfr, label='input', dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)
    srwlpy.SetRepresElecField(wfr, 'f')

    t0 = time(); print('Resizing in frequency domain: ', end='')
    srwlpy.ResizeElecField(wfr, 'f', [0, 1., z_scaling]); print('done in', round(time() - t0, 3), 's')

    # diagnostics - input
    diagnostics, diagnostics_names = diagnose_input(wfr)

    # HHLM
    label = 'before HHLM1'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_CRL0_before_HHLM1(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before HHLM2'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_HHLM1_before_HHLM2(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before HHLM3'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_HHLM2_before_HHLM3(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before HHLM4'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_HHLM3_before_HHLM4(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'output'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_HHLM4_output(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    print('\n\neverything lasted: {}s'.format(round(time()-time_stamp,2)))

    # diagnostics - output
    diagnostics, diagnostics_names = diagnose_output(wfr, diagnostics=diagnostics, diagnostics_names=diagnostics_names)

    # I/O
    dir_diagnostics = dir_plot+'diagnostics/'; mkdir(dir_diagnostics)
    save_diagnostics(dir_diagnostics, diagnostics, diagnostics_names)

    print('beam diagnostics written to {}'.format(dir_diagnostics))


#### HRM
def main_HRM(varParam_name, z_scaling=10.0, if_close=0, dir_plot=None, if_log=1, if_slice=1, i_start=0):
    time_stamp=time()

    # initialization
    with open(varParam_name, 'rb') as f:
        varParam = pickle.load(f)
    v = srwl_bl.srwl_uti_parse_options(varParam, use_sys_argv=True)
    v.si = True
    srwl_bl.SRWLBeamline(_name=v.name).calc_all(v)

    # incident beam
    wfr = v.w_res; i_plot = 1
    plot_wfr_diagnostic(wfr, label='input', dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)
    srwlpy.SetRepresElecField(wfr, 'f')

    t0 = time(); print('Resizing in frequency domain: ', end='')
    srwlpy.ResizeElecField(wfr, 'f', [0, 1., z_scaling]); print('done in', round(time() - t0, 3), 's')

    # diagnostics - input
    diagnostics, diagnostics_names = diagnose_input(wfr)

    # CC1
    label = 'before C1'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_CRL0_before_C1(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before C2'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_C1_before_C2(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before CRL1'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_C2_before_CRL1(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    # Telescope
    label_append = ', open'
    label = 'focal plane'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_CRL1_focus(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    if if_close == 1:
        label_append = ', closed'
        label = 'focal plane'+label_append; print('Propagating to slit: ', end=''); t0 = time()
        bl = set_optics_slit(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before CRL2'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_slit_before_CRL2(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before C3'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_CRL2_before_C3(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    # CC2
    label = 'before C4'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_C3_before_C4(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'output'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_C4_output(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    print('\n\neverything lasted: {}s'.format(round(time()-time_stamp,2)))

    # diagnostics - output
    diagnostics, diagnostics_names = diagnose_output(wfr, diagnostics=diagnostics, diagnostics_names=diagnostics_names)

    # I/O
    if if_close == 1:
        dir_diagnostics = dir_plot+'diagnostics_closed/'; mkdir(dir_diagnostics)
    else:
        dir_diagnostics = dir_plot+'diagnostics_open/'; mkdir(dir_diagnostics)
    save_diagnostics(dir_diagnostics, diagnostics, diagnostics_names)

    print('\nbeam diagnostics written to {}'.format(dir_diagnostics))


#### HHLM+HRM
def main_HHLM_HRM(varParam_name, z_scaling=10.0, if_close=0, dir_plot=None, if_log=1, if_slice=1, i_start=0):
    time_stamp=time()

    # initialization
    with open(varParam_name, 'rb') as f:
        varParam = pickle.load(f)
    v = srwl_bl.srwl_uti_parse_options(varParam, use_sys_argv=True)
    v.si = True
    srwl_bl.SRWLBeamline(_name=v.name).calc_all(v)

    # incident beam
    wfr = v.w_res; i_plot = 1
    plot_wfr_diagnostic(wfr, label='input', dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)
    srwlpy.SetRepresElecField(wfr, 'f')

    t0 = time(); print('Resizing in frequency domain: ', end='')
    srwlpy.ResizeElecField(wfr, 'f', [0, 1., z_scaling]); print('done in', round(time() - t0, 3), 's')

    # diagnostics - input
    diagnostics, diagnostics_names = diagnose_input(wfr)

    # HHLM
    label = 'before HHLM1'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_CRL0_before_HHLM1(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before HHLM2'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_HHLM1_before_HHLM2(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before HHLM3'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_HHLM2_before_HHLM3(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before HHLM4'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_HHLM3_before_HHLM4(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before C1'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_CRL0_before_C1(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    # CC1
    label = 'before C2'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_C1_before_C2(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before CRL1'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_C2_before_CRL1(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    # Telescope
    label_append = ', open'
    label = 'focal plane'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_CRL1_focus(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    if if_close == 1:
        label_append = ', closed'
        label = 'focal plane'+label_append; print('Propagating to slit: ', end=''); t0 = time()
        bl = set_optics_slit(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before CRL2'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_slit_before_CRL2(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'before C3'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_CRL2_before_C3(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    # CC2
    label = 'before C4'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_C3_before_C4(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    label = 'output'+label_append; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    bl = set_optics_C4_output(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    print('\n\neverything lasted: {}s'.format(round(time()-time_stamp,2)))

    # diagnostics - output
    diagnostics, diagnostics_names = diagnose_output(wfr, diagnostics=diagnostics, diagnostics_names=diagnostics_names)

    # I/O
    if if_close == 1:
        dir_diagnostics = dir_plot+'diagnostics_closed/'; mkdir(dir_diagnostics)
    else:
        dir_diagnostics = dir_plot+'diagnostics_open/'; mkdir(dir_diagnostics)
    save_diagnostics(dir_diagnostics, diagnostics, diagnostics_names)

    print('\nbeam diagnostics written to {}'.format(dir_diagnostics))
