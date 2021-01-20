#!/usr/bin/env python
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
import math
import srwl_uti_smp


def set_optics(v=None):
    el = []
    pp = []
    names = ['MR1L3', 'MR1L3_MR2L3', 'MR2L3', 'MR2L3_Before_HHLM', 'C1', 'After_C1_C2', 'C2', 'After_C2_C3', 'C3', 'After_C3_C4', 'C4', 'After_C4', 'After_C4_HRM1', 'HRM1', 'After_HRM1_HRM2', 'HRM2', 'After_HRM2', 'After_HRM2_HRM_F1', 'HRM_F1', 'After_HRM_F1', 'After_HRM_F1_Before_WDS', 'WDS', 'After_WDS', 'After_WDS_HRM_F2', 'HRM_F2', 'After_HRM_F2', 'After_HRM_F2_Before_HRM3', 'Before_HRM3', 'HRM3', 'After_HRM3', 'After_HRM3_HRM4', 'HRM4', 'After_HRM4', 'After_HRM4_After_HRM', 'After_HRM']
    for el_name in names:
        if el_name == 'MR1L3':
            # MR1L3: mirror 157.05m
            mirror_file = v.op_MR1L3_hfn
            assert os.path.isfile(mirror_file), \
                'Missing input file {}, required by MR1L3 beamline element'.format(mirror_file)
            el.append(srwlib.srwl_opt_setup_surf_height_1d(
                srwlib.srwl_uti_read_data_cols(mirror_file, "\t", 0, 1),
                _dim=v.op_MR1L3_dim,
                _ang=abs(v.op_MR1L3_ang),
                _amp_coef=v.op_MR1L3_amp_coef,
                _size_x=v.op_MR1L3_size_x,
                _size_y=v.op_MR1L3_size_y,
            ))
            pp.append(v.op_MR1L3_pp)
        elif el_name == 'MR1L3_MR2L3':
            # MR1L3_MR2L3: drift 157.05m
            el.append(srwlib.SRWLOptD(
                _L=v.op_MR1L3_MR2L3_L,
            ))
            pp.append(v.op_MR1L3_MR2L3_pp)
        elif el_name == 'MR2L3':
            # MR2L3: mirror 271.05m
            mirror_file = v.op_MR2L3_hfn
            assert os.path.isfile(mirror_file), \
                'Missing input file {}, required by MR2L3 beamline element'.format(mirror_file)
            el.append(srwlib.srwl_opt_setup_surf_height_1d(
                srwlib.srwl_uti_read_data_cols(mirror_file, "\t", 0, 1),
                _dim=v.op_MR2L3_dim,
                _ang=abs(v.op_MR2L3_ang),
                _amp_coef=v.op_MR2L3_amp_coef,
                _size_x=v.op_MR2L3_size_x,
                _size_y=v.op_MR2L3_size_y,
            ))
            pp.append(v.op_MR2L3_pp)
        elif el_name == 'MR2L3_Before_HHLM':
            # MR2L3_Before_HHLM: drift 271.05m
            el.append(srwlib.SRWLOptD(
                _L=v.op_MR2L3_Before_HHLM_L,
            ))
            pp.append(v.op_MR2L3_Before_HHLM_pp)
        elif el_name == 'C1':
            # C1: crystal 279.05m
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

        elif el_name == 'After_C1_C2':
            # After_C1_C2: drift 279.05m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_C1_C2_L,
            ))
            pp.append(v.op_After_C1_C2_pp)
        elif el_name == 'C2':
            # C2: crystal 279.09m
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

        elif el_name == 'After_C2_C3':
            # After_C2_C3: drift 279.09m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_C2_C3_L,
            ))
            pp.append(v.op_After_C2_C3_pp)
        elif el_name == 'C3':
            # C3: crystal 279.34m
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

        elif el_name == 'After_C3_C4':
            # After_C3_C4: drift 279.34m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_C3_C4_L,
            ))
            pp.append(v.op_After_C3_C4_pp)
        elif el_name == 'C4':
            # C4: crystal 279.55m
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

        elif el_name == 'After_C4':
            # After_C4: watch 279.55m
            pass
        elif el_name == 'After_C4_HRM1':
            # After_C4_HRM1: drift 279.55m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_C4_HRM1_L,
            ))
            pp.append(v.op_After_C4_HRM1_pp)
        elif el_name == 'HRM1':
            # HRM1: crystal 289.47m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_HRM1_d_sp,
                _psi0r=v.op_HRM1_psi0r,
                _psi0i=v.op_HRM1_psi0i,
                _psi_hr=v.op_HRM1_psiHr,
                _psi_hi=v.op_HRM1_psiHi,
                _psi_hbr=v.op_HRM1_psiHBr,
                _psi_hbi=v.op_HRM1_psiHBi,
                _tc=v.op_HRM1_tc,
                _ang_as=v.op_HRM1_ang_as,
                _nvx=v.op_HRM1_nvx,
                _nvy=v.op_HRM1_nvy,
                _nvz=v.op_HRM1_nvz,
                _tvx=v.op_HRM1_tvx,
                _tvy=v.op_HRM1_tvy,
                _uc=v.op_HRM1_uc,
                _e_avg=v.op_HRM1_energy,
                _ang_roll=v.op_HRM1_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_HRM1_pp)

        elif el_name == 'After_HRM1_HRM2':
            # After_HRM1_HRM2: drift 289.47m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_HRM1_HRM2_L,
            ))
            pp.append(v.op_After_HRM1_HRM2_pp)
        elif el_name == 'HRM2':
            # HRM2: crystal 289.55m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_HRM2_d_sp,
                _psi0r=v.op_HRM2_psi0r,
                _psi0i=v.op_HRM2_psi0i,
                _psi_hr=v.op_HRM2_psiHr,
                _psi_hi=v.op_HRM2_psiHi,
                _psi_hbr=v.op_HRM2_psiHBr,
                _psi_hbi=v.op_HRM2_psiHBi,
                _tc=v.op_HRM2_tc,
                _ang_as=v.op_HRM2_ang_as,
                _nvx=v.op_HRM2_nvx,
                _nvy=v.op_HRM2_nvy,
                _nvz=v.op_HRM2_nvz,
                _tvx=v.op_HRM2_tvx,
                _tvy=v.op_HRM2_tvy,
                _uc=v.op_HRM2_uc,
                _e_avg=v.op_HRM2_energy,
                _ang_roll=v.op_HRM2_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_HRM2_pp)

        elif el_name == 'After_HRM2':
            # After_HRM2: watch 289.55m
            pass
        elif el_name == 'After_HRM2_HRM_F1':
            # After_HRM2_HRM_F1: drift 289.55m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_HRM2_HRM_F1_L,
            ))
            pp.append(v.op_After_HRM2_HRM_F1_pp)
        elif el_name == 'HRM_F1':
            # HRM_F1: crl 299.55m
            el.append(srwlib.srwl_opt_setup_CRL(
                _foc_plane=v.op_HRM_F1_foc_plane,
                _delta=v.op_HRM_F1_delta,
                _atten_len=v.op_HRM_F1_atten_len,
                _shape=v.op_HRM_F1_shape,
                _apert_h=v.op_HRM_F1_apert_h,
                _apert_v=v.op_HRM_F1_apert_v,
                _r_min=v.op_HRM_F1_r_min,
                _n=v.op_HRM_F1_n,
                _wall_thick=v.op_HRM_F1_wall_thick,
                _xc=v.op_HRM_F1_x,
                _yc=v.op_HRM_F1_y,
            ))
            pp.append(v.op_HRM_F1_pp)
        elif el_name == 'After_HRM_F1':
            # After_HRM_F1: watch 299.55m
            pass
        elif el_name == 'After_HRM_F1_Before_WDS':
            # After_HRM_F1_Before_WDS: drift 299.55m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_HRM_F1_Before_WDS_L,
            ))
            pp.append(v.op_After_HRM_F1_Before_WDS_pp)
        elif el_name == 'WDS':
            # WDS: aperture 309.55m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_WDS_shape,
                _ap_or_ob='a',
                _Dx=v.op_WDS_Dx,
                _Dy=v.op_WDS_Dy,
                _x=v.op_WDS_x,
                _y=v.op_WDS_y,
            ))
            pp.append(v.op_WDS_pp)
        elif el_name == 'After_WDS':
            # After_WDS: watch 309.55m
            pass
        elif el_name == 'After_WDS_HRM_F2':
            # After_WDS_HRM_F2: drift 309.55m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_WDS_HRM_F2_L,
            ))
            pp.append(v.op_After_WDS_HRM_F2_pp)
        elif el_name == 'HRM_F2':
            # HRM_F2: crl 319.55m
            el.append(srwlib.srwl_opt_setup_CRL(
                _foc_plane=v.op_HRM_F2_foc_plane,
                _delta=v.op_HRM_F2_delta,
                _atten_len=v.op_HRM_F2_atten_len,
                _shape=v.op_HRM_F2_shape,
                _apert_h=v.op_HRM_F2_apert_h,
                _apert_v=v.op_HRM_F2_apert_v,
                _r_min=v.op_HRM_F2_r_min,
                _n=v.op_HRM_F2_n,
                _wall_thick=v.op_HRM_F2_wall_thick,
                _xc=v.op_HRM_F2_x,
                _yc=v.op_HRM_F2_y,
            ))
            pp.append(v.op_HRM_F2_pp)
        elif el_name == 'After_HRM_F2':
            # After_HRM_F2: watch 319.55m
            pass
        elif el_name == 'After_HRM_F2_Before_HRM3':
            # After_HRM_F2_Before_HRM3: drift 319.55m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_HRM_F2_Before_HRM3_L,
            ))
            pp.append(v.op_After_HRM_F2_Before_HRM3_pp)
        elif el_name == 'Before_HRM3':
            # Before_HRM3: watch 329.55m
            pass
        elif el_name == 'HRM3':
            # HRM3: crystal 329.55m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_HRM3_d_sp,
                _psi0r=v.op_HRM3_psi0r,
                _psi0i=v.op_HRM3_psi0i,
                _psi_hr=v.op_HRM3_psiHr,
                _psi_hi=v.op_HRM3_psiHi,
                _psi_hbr=v.op_HRM3_psiHBr,
                _psi_hbi=v.op_HRM3_psiHBi,
                _tc=v.op_HRM3_tc,
                _ang_as=v.op_HRM3_ang_as,
                _nvx=v.op_HRM3_nvx,
                _nvy=v.op_HRM3_nvy,
                _nvz=v.op_HRM3_nvz,
                _tvx=v.op_HRM3_tvx,
                _tvy=v.op_HRM3_tvy,
                _uc=v.op_HRM3_uc,
                _e_avg=v.op_HRM3_energy,
                _ang_roll=v.op_HRM3_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_HRM3_pp)

        elif el_name == 'After_HRM3':
            # After_HRM3: watch 329.55m
            pass
        elif el_name == 'After_HRM3_HRM4':
            # After_HRM3_HRM4: drift 329.55m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_HRM3_HRM4_L,
            ))
            pp.append(v.op_After_HRM3_HRM4_pp)
        elif el_name == 'HRM4':
            # HRM4: crystal 329.63m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_HRM4_d_sp,
                _psi0r=v.op_HRM4_psi0r,
                _psi0i=v.op_HRM4_psi0i,
                _psi_hr=v.op_HRM4_psiHr,
                _psi_hi=v.op_HRM4_psiHi,
                _psi_hbr=v.op_HRM4_psiHBr,
                _psi_hbi=v.op_HRM4_psiHBi,
                _tc=v.op_HRM4_tc,
                _ang_as=v.op_HRM4_ang_as,
                _nvx=v.op_HRM4_nvx,
                _nvy=v.op_HRM4_nvy,
                _nvz=v.op_HRM4_nvz,
                _tvx=v.op_HRM4_tvx,
                _tvy=v.op_HRM4_tvy,
                _uc=v.op_HRM4_uc,
                _e_avg=v.op_HRM4_energy,
                _ang_roll=v.op_HRM4_diffractionAngle
            )
            el.append(crystal)
            pp.append(v.op_HRM4_pp)

        elif el_name == 'After_HRM4':
            # After_HRM4: watch 329.63m
            pass
        elif el_name == 'After_HRM4_After_HRM':
            # After_HRM4_After_HRM: drift 329.63m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_HRM4_After_HRM_L,
            ))
            pp.append(v.op_After_HRM4_After_HRM_pp)
        elif el_name == 'After_HRM':
            # After_HRM: watch 339.63m
            pass
    pp.append(v.op_fin_pp)
    return srwlib.SRWLOptC(el, pp)


varParam = srwl_bl.srwl_uti_ext_options([
    ['name', 's', 'Gaussian Beam 9481 eV #01', 'simulation name'],

#---Data Folder
    ['fdir', 's', '', 'folder (directory) name for reading-in input and saving output data files'],


    ['gbm_x', 'f', 0.0, 'average horizontal coordinates of waist [m]'],
    ['gbm_y', 'f', 0.0, 'average vertical coordinates of waist [m]'],
    ['gbm_z', 'f', 0.0, 'average longitudinal coordinate of waist [m]'],
    ['gbm_xp', 'f', 0.0, 'average horizontal angle at waist [rad]'],
    ['gbm_yp', 'f', 0.0, 'average verical angle at waist [rad]'],
    ['gbm_ave', 'f', 9481.0, 'average photon energy [eV]'],
    ['gbm_pen', 'f', 0.001, 'energy per pulse [J]'],
    ['gbm_rep', 'f', 1, 'rep. rate [Hz]'],
    ['gbm_pol', 'f', 2, 'polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left'],
    ['gbm_sx', 'f', 1.4792385e-05, 'rms beam size vs horizontal position [m] at waist (for intensity)'],
    ['gbm_sy', 'f', 1.528112e-05, 'rms beam size vs vertical position [m] at waist (for intensity)'],
    ['gbm_st', 'f', 1e-13, 'rms pulse duration [s] (for intensity)'],
    ['gbm_mx', 'f', 0, 'transverse Gauss-Hermite mode order in horizontal direction'],
    ['gbm_my', 'f', 0, 'transverse Gauss-Hermite mode order in vertical direction'],
    ['gbm_ca', 's', 'c', 'treat _sigX, _sigY as sizes in [m] in coordinate representation (_presCA="c") or as angular divergences in [rad] in angular representation (_presCA="a")'],
    ['gbm_ft', 's', 't', 'treat _sigT as pulse duration in [s] in time domain/representation (_presFT="t") or as bandwidth in [eV] in frequency domain/representation (_presFT="f")'],

#---Calculation Types
    # Electron Trajectory
    ['tr', '', '', 'calculate electron trajectory', 'store_true'],
    ['tr_cti', 'f', 0.0, 'initial time moment (c*t) for electron trajectory calculation [m]'],
    ['tr_ctf', 'f', 0.0, 'final time moment (c*t) for electron trajectory calculation [m]'],
    ['tr_np', 'f', 10000, 'number of points for trajectory calculation'],
    ['tr_mag', 'i', 1, 'magnetic field to be used for trajectory calculation: 1- approximate, 2- accurate'],
    ['tr_fn', 's', 'res_trj.dat', 'file name for saving calculated trajectory data'],
    ['tr_pl', 's', '', 'plot the resulting trajectiry in graph(s): ""- dont plot, otherwise the string should list the trajectory components to plot'],

    #Single-Electron Spectrum vs Photon Energy
    ['ss', '', '', 'calculate single-e spectrum vs photon energy', 'store_true'],
    ['ss_ei', 'f', 100.0, 'initial photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ef', 'f', 20000.0, 'final photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ne', 'i', 10000, 'number of points vs photon energy for single-e spectrum vs photon energy calculation'],
    ['ss_x', 'f', 0.0, 'horizontal position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_y', 'f', 0.0, 'vertical position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_meth', 'i', 1, 'method to use for single-e spectrum vs photon energy calculation: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['ss_prec', 'f', 0.01, 'relative precision for single-e spectrum vs photon energy calculation (nominal value is 0.01)'],
    ['ss_pol', 'i', 6, 'polarization component to extract after spectrum vs photon energy calculation: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['ss_mag', 'i', 1, 'magnetic field to be used for single-e spectrum vs photon energy calculation: 1- approximate, 2- accurate'],
    ['ss_ft', 's', 'f', 'presentation/domain: "f"- frequency (photon energy), "t"- time'],
    ['ss_u', 'i', 1, 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['ss_fn', 's', 'res_spec_se.dat', 'file name for saving calculated single-e spectrum vs photon energy'],
    ['ss_pl', 's', '', 'plot the resulting single-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],

    #Multi-Electron Spectrum vs Photon Energy (taking into account e-beam emittance, energy spread and collection aperture size)
    ['sm', '', '', 'calculate multi-e spectrum vs photon energy', 'store_true'],
    ['sm_ei', 'f', 100.0, 'initial photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ef', 'f', 20000.0, 'final photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ne', 'i', 10000, 'number of points vs photon energy for multi-e spectrum vs photon energy calculation'],
    ['sm_x', 'f', 0.0, 'horizontal center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_rx', 'f', 0.001, 'range of horizontal position / horizontal aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_nx', 'i', 1, 'number of points vs horizontal position for multi-e spectrum vs photon energy calculation'],
    ['sm_y', 'f', 0.0, 'vertical center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ry', 'f', 0.001, 'range of vertical position / vertical aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ny', 'i', 1, 'number of points vs vertical position for multi-e spectrum vs photon energy calculation'],
    ['sm_mag', 'i', 1, 'magnetic field to be used for calculation of multi-e spectrum spectrum or intensity distribution: 1- approximate, 2- accurate'],
    ['sm_hi', 'i', 1, 'initial UR spectral harmonic to be taken into account for multi-e spectrum vs photon energy calculation'],
    ['sm_hf', 'i', 15, 'final UR spectral harmonic to be taken into account for multi-e spectrum vs photon energy calculation'],
    ['sm_prl', 'f', 1.0, 'longitudinal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_pra', 'f', 1.0, 'azimuthal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_meth', 'i', -1, 'method to use for spectrum vs photon energy calculation in case of arbitrary input magnetic field: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler", -1- dont use this accurate integration method (rather use approximate if possible)'],
    ['sm_prec', 'f', 0.01, 'relative precision for spectrum vs photon energy calculation in case of arbitrary input magnetic field (nominal value is 0.01)'],
    ['sm_nm', 'i', 1, 'number of macro-electrons for calculation of spectrum in case of arbitrary input magnetic field'],
    ['sm_na', 'i', 5, 'number of macro-electrons to average on each node at parallel (MPI-based) calculation of spectrum in case of arbitrary input magnetic field'],
    ['sm_ns', 'i', 5, 'saving periodicity (in terms of macro-electrons) for intermediate intensity at calculation of multi-electron spectrum in case of arbitrary input magnetic field'],
    ['sm_type', 'i', 1, 'calculate flux (=1) or flux per unit surface (=2)'],
    ['sm_pol', 'i', 6, 'polarization component to extract after calculation of multi-e flux or intensity: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['sm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['sm_fn', 's', 'res_spec_me.dat', 'file name for saving calculated milti-e spectrum vs photon energy'],
    ['sm_pl', 's', '', 'plot the resulting spectrum-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],
    #to add options for the multi-e calculation from "accurate" magnetic field

    #Power Density Distribution vs horizontal and vertical position
    ['pw', '', '', 'calculate SR power density distribution', 'store_true'],
    ['pw_x', 'f', 0.0, 'central horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_rx', 'f', 0.015, 'range of horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_nx', 'i', 100, 'number of points vs horizontal position for calculation of power density distribution'],
    ['pw_y', 'f', 0.0, 'central vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ry', 'f', 0.015, 'range of vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ny', 'i', 100, 'number of points vs vertical position for calculation of power density distribution'],
    ['pw_pr', 'f', 1.0, 'precision factor for calculation of power density distribution'],
    ['pw_meth', 'i', 1, 'power density computation method (1- "near field", 2- "far field")'],
    ['pw_zst', 'f', 0., 'initial longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_zfi', 'f', 0., 'final longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_mag', 'i', 1, 'magnetic field to be used for power density calculation: 1- approximate, 2- accurate'],
    ['pw_fn', 's', 'res_pow.dat', 'file name for saving calculated power density distribution'],
    ['pw_pl', 's', '', 'plot the resulting power density distribution in a graph: ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    #Single-Electron Intensity distribution vs horizontal and vertical position
    ['si', '', '', 'calculate single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position', 'store_true'],
    #Single-Electron Wavefront Propagation
    ['ws', '', '', 'calculate single-electron (/ fully coherent) wavefront propagation', 'store_true'],
    #Multi-Electron (partially-coherent) Wavefront Propagation
    ['wm', '', '', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],

    ['w_e', 'f', 9481.0, 'photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ef', 'f', -1.0, 'final photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ne', 'i', 1, 'number of points vs photon energy for calculation of intensity distribution'],
    ['w_x', 'f', 0.0, 'central horizontal position [m] for calculation of intensity distribution'],
    ['w_rx', 'f', 0.00265, 'range of horizontal position [m] for calculation of intensity distribution'],
    ['w_nx', 'i', 100, 'number of points vs horizontal position for calculation of intensity distribution'],
    ['w_y', 'f', 0.0, 'central vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ry', 'f', 0.00265, 'range of vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ny', 'i', 100, 'number of points vs vertical position for calculation of intensity distribution'],
    ['w_smpf', 'f', 0.5, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_meth', 'i', 2, 'method to use for calculation of intensity distribution vs horizontal and vertical position: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['w_prec', 'f', 0.01, 'relative precision for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_u', 'i', 1, 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['si_pol', 'i', 6, 'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['si_type', 'i', 0, 'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],
    ['w_mag', 'i', 1, 'magnetic field to be used for calculation of intensity distribution vs horizontal and vertical position: 1- approximate, 2- accurate'],

    ['si_fn', 's', 'res_int_se.dat', 'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],
    ['si_pl', 's', '', 'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
    ['ws_fni', 's', 'res_int_pr_se.dat', 'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],
    ['ws_pl', 's', '', 'plot the resulting intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    ['wm_nm', 'i', 1000, 'number of macro-electrons (coherent wavefronts) for calculation of multi-electron wavefront propagation'],
    ['wm_na', 'i', 5, 'number of macro-electrons (coherent wavefronts) to average on each node for parallel (MPI-based) calculation of multi-electron wavefront propagation'],
    ['wm_ns', 'i', 5, 'saving periodicity (in terms of macro-electrons / coherent wavefronts) for intermediate intensity at multi-electron wavefront propagation calculation'],
    ['wm_ch', 'i', 0, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y; 40- intensity(s0), mutual intensity cuts and degree of coherence vs X & Y'],
    ['wm_ap', 'i', 0, 'switch specifying representation of the resulting Stokes parameters: coordinate (0) or angular (1)'],
    ['wm_x0', 'f', 0.0, 'horizontal center position for mutual intensity cut calculation'],
    ['wm_y0', 'f', 0.0, 'vertical center position for mutual intensity cut calculation'],
    ['wm_ei', 'i', 0, 'integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from w_e, w_ef'],
    ['wm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['wm_am', 'i', 0, 'multi-electron integration approximation method: 0- no approximation (use the standard 5D integration method), 1- integrate numerically only over e-beam energy spread and use convolution to treat transverse emittance'],
    ['wm_fni', 's', 'res_int_pr_me.dat', 'file name for saving propagated multi-e intensity distribution vs horizontal and vertical position'],
    ['wm_fbk', '', '', 'create backup file(s) with propagated multi-e intensity distribution vs horizontal and vertical position and other radiation characteristics', 'store_true'],

    #to add options
    ['op_r', 'f', 157.05, 'longitudinal position of the first optical element [m]'],
    # Former appParam:
    ['rs_type', 's', 'g', 'source type, (u) idealized undulator, (t), tabulated undulator, (m) multipole, (g) gaussian beam'],

#---Beamline optics:
    # MR1L3: mirror
    ['op_MR1L3_hfn', 's', 'mirror_1d.dat', 'heightProfileFile'],
    ['op_MR1L3_dim', 's', 'x', 'orientation'],
    ['op_MR1L3_ang', 'f', 0.00265, 'grazingAngle'],
    ['op_MR1L3_amp_coef', 'f', 0.001, 'heightAmplification'],
    ['op_MR1L3_size_x', 'f', 0.00265, 'horizontalTransverseSize'],
    ['op_MR1L3_size_y', 'f', 0.01, 'verticalTransverseSize'],

    # MR1L3_MR2L3: drift
    ['op_MR1L3_MR2L3_L', 'f', 114.0, 'length'],

    # MR2L3: mirror
    ['op_MR2L3_hfn', 's', 'mirror_1d.dat', 'heightProfileFile'],
    ['op_MR2L3_dim', 's', 'x', 'orientation'],
    ['op_MR2L3_ang', 'f', -0.00265, 'grazingAngle'],
    ['op_MR2L3_amp_coef', 'f', 0.001, 'heightAmplification'],
    ['op_MR2L3_size_x', 'f', 0.00265, 'horizontalTransverseSize'],
    ['op_MR2L3_size_y', 'f', 0.01, 'verticalTransverseSize'],

    # MR2L3_Before_HHLM: drift
    ['op_MR2L3_Before_HHLM_L', 'f', 8.0, 'length'],

    # C1: crystal
    ['op_C1_hfn', 's', '', 'heightProfileFile'],
    ['op_C1_dim', 's', 'x', 'orientation'],
    ['op_C1_d_sp', 'f', 1.9201374688016222, 'dSpacing'],
    ['op_C1_psi0r', 'f', -1.0873035035585694e-05, 'psi0r'],
    ['op_C1_psi0i', 'f', 1.8438837339536554e-07, 'psi0i'],
    ['op_C1_psiHr', 'f', -6.610167427198717e-06, 'psiHr'],
    ['op_C1_psiHi', 'f', 1.7824173540780476e-07, 'psiHi'],
    ['op_C1_psiHBr', 'f', -6.610167427198717e-06, 'psiHBr'],
    ['op_C1_psiHBi', 'f', 1.7824173540780476e-07, 'psiHBi'],
    ['op_C1_tc', 'f', 0.01, 'crystalThickness'],
    ['op_C1_uc', 'f', 1, 'useCase'],
    ['op_C1_ang_as', 'f', 0.296706, 'asymmetryAngle'],
    ['op_C1_nvx', 'f', -0.9987059480792639, 'nvx'],
    ['op_C1_nvy', 'f', 6.786103757443706e-09, 'nvy'],
    ['op_C1_nvz', 'f', -0.05085694909349798, 'nvz'],
    ['op_C1_tvx', 'f', -0.05085694909349798, 'tvx'],
    ['op_C1_tvy', 'f', 3.455677159020173e-10, 'tvy'],
    ['op_C1_ang', 'f', 0.05087889763248938, 'grazingAngle'],
    ['op_C1_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_C1_energy', 'f', 9481.0, 'energy'],
    ['op_C1_diffractionAngle', 'f', 1.57079632, 'diffractionAngle'],

    # After_C1_C2: drift
    ['op_After_C1_C2_L', 'f', 0.03999999999996362, 'length'],

    # C2: crystal
    ['op_C2_hfn', 's', '', 'heightProfileFile'],
    ['op_C2_dim', 's', 'x', 'orientation'],
    ['op_C2_d_sp', 'f', 0.9600687344008111, 'dSpacing'],
    ['op_C2_psi0r', 'f', -1.0873035035585694e-05, 'psi0r'],
    ['op_C2_psi0i', 'f', 1.8438837339536554e-07, 'psi0i'],
    ['op_C2_psiHr', 'f', -4.181686438547451e-06, 'psiHr'],
    ['op_C2_psiHi', 'f', 1.6100412693351052e-07, 'psiHi'],
    ['op_C2_psiHBr', 'f', -4.181686438547451e-06, 'psiHBr'],
    ['op_C2_psiHBi', 'f', 1.6100412693351052e-07, 'psiHBi'],
    ['op_C2_tc', 'f', 0.01, 'crystalThickness'],
    ['op_C2_uc', 'f', 1, 'useCase'],
    ['op_C2_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_C2_nvx', 'f', 0.7322282430733594, 'nvx'],
    ['op_C2_nvy', 'f', 4.975415277322606e-09, 'nvy'],
    ['op_C2_nvz', 'f', -0.6810593219725439, 'nvz'],
    ['op_C2_tvx', 'f', 0.6810593219725439, 'tvx'],
    ['op_C2_tvy', 'f', 4.627727743855522e-09, 'tvy'],
    ['op_C2_ang', 'f', 0.7492083731847909, 'grazingAngle'],
    ['op_C2_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_C2_energy', 'f', 9481.0, 'energy'],
    ['op_C2_diffractionAngle', 'f', -1.57079632, 'diffractionAngle'],

    # After_C2_C3: drift
    ['op_After_C2_C3_L', 'f', 0.25, 'length'],

    # C3: crystal
    ['op_C3_hfn', 's', '', 'heightProfileFile'],
    ['op_C3_dim', 's', 'x', 'orientation'],
    ['op_C3_d_sp', 'f', 0.9600687344008111, 'dSpacing'],
    ['op_C3_psi0r', 'f', -1.0873035035585694e-05, 'psi0r'],
    ['op_C3_psi0i', 'f', 1.8438837339536554e-07, 'psi0i'],
    ['op_C3_psiHr', 'f', -4.181686438547451e-06, 'psiHr'],
    ['op_C3_psiHi', 'f', 1.6100412693351052e-07, 'psiHi'],
    ['op_C3_psiHBr', 'f', -4.181686438547451e-06, 'psiHBr'],
    ['op_C3_psiHBi', 'f', 1.6100412693351052e-07, 'psiHBi'],
    ['op_C3_tc', 'f', 0.01, 'crystalThickness'],
    ['op_C3_uc', 'f', 1, 'useCase'],
    ['op_C3_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_C3_nvx', 'f', -0.7322282430733594, 'nvx'],
    ['op_C3_nvy', 'f', 4.975415277322606e-09, 'nvy'],
    ['op_C3_nvz', 'f', -0.6810593219725439, 'nvz'],
    ['op_C3_tvx', 'f', -0.6810593219725439, 'tvx'],
    ['op_C3_tvy', 'f', 4.627727743855522e-09, 'tvy'],
    ['op_C3_ang', 'f', 0.7492083731847909, 'grazingAngle'],
    ['op_C3_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_C3_energy', 'f', 9481.0, 'energy'],
    ['op_C3_diffractionAngle', 'f', 1.57079632, 'diffractionAngle'],

    # After_C3_C4: drift
    ['op_After_C3_C4_L', 'f', 0.21000000000003638, 'length'],

    # C4: crystal
    ['op_C4_hfn', 's', '', 'heightProfileFile'],
    ['op_C4_dim', 's', 'x', 'orientation'],
    ['op_C4_d_sp', 'f', 1.9201374688016222, 'dSpacing'],
    ['op_C4_psi0r', 'f', -1.0873035035585694e-05, 'psi0r'],
    ['op_C4_psi0i', 'f', 1.8438837339536554e-07, 'psi0i'],
    ['op_C4_psiHr', 'f', -6.610167427198717e-06, 'psiHr'],
    ['op_C4_psiHi', 'f', 1.7824173540780476e-07, 'psiHi'],
    ['op_C4_psiHBr', 'f', -6.610167427198717e-06, 'psiHBr'],
    ['op_C4_psiHBi', 'f', 1.7824173540780476e-07, 'psiHBi'],
    ['op_C4_tc', 'f', 0.01, 'crystalThickness'],
    ['op_C4_uc', 'f', 1, 'useCase'],
    ['op_C4_ang_as', 'f', -0.296706, 'asymmetryAngle'],
    ['op_C4_nvx', 'f', 0.7995857954948998, 'nvx'],
    ['op_C4_nvy', 'f', 5.43310288843489e-09, 'nvy'],
    ['op_C4_nvz', 'f', -0.6005518758964858, 'nvz'],
    ['op_C4_tvx', 'f', 0.6005518758964858, 'tvx'],
    ['op_C4_tvy', 'f', 4.0806879636583054e-09, 'tvy'],
    ['op_C4_ang', 'f', 0.6441911322682903, 'grazingAngle'],
    ['op_C4_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_C4_energy', 'f', 9481.0, 'energy'],
    ['op_C4_diffractionAngle', 'f', -1.57079632, 'diffractionAngle'],

    # After_C4_HRM1: drift
    ['op_After_C4_HRM1_L', 'f', 9.920000000000016, 'length'],

    # HRM1: crystal
    ['op_HRM1_hfn', 's', '', 'heightProfileFile'],
    ['op_HRM1_dim', 's', 'x', 'orientation'],
    ['op_HRM1_d_sp', 'f', 0.9600687344008111, 'dSpacing'],
    ['op_HRM1_psi0r', 'f', -1.0873035035585694e-05, 'psi0r'],
    ['op_HRM1_psi0i', 'f', 1.8438837339536554e-07, 'psi0i'],
    ['op_HRM1_psiHr', 'f', -4.181686438547451e-06, 'psiHr'],
    ['op_HRM1_psiHi', 'f', 1.6100412693351052e-07, 'psiHi'],
    ['op_HRM1_psiHBr', 'f', -4.181686438547451e-06, 'psiHBr'],
    ['op_HRM1_psiHBi', 'f', 1.6100412693351052e-07, 'psiHBi'],
    ['op_HRM1_tc', 'f', 0.01, 'crystalThickness'],
    ['op_HRM1_uc', 'f', 1, 'useCase'],
    ['op_HRM1_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_HRM1_nvx', 'f', -0.7322282430733594, 'nvx'],
    ['op_HRM1_nvy', 'f', 4.975415277322606e-09, 'nvy'],
    ['op_HRM1_nvz', 'f', -0.6810593219725439, 'nvz'],
    ['op_HRM1_tvx', 'f', -0.6810593219725439, 'tvx'],
    ['op_HRM1_tvy', 'f', 4.627727743855522e-09, 'tvy'],
    ['op_HRM1_ang', 'f', 0.7492083731847909, 'grazingAngle'],
    ['op_HRM1_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_HRM1_energy', 'f', 9481.0, 'energy'],
    ['op_HRM1_diffractionAngle', 'f', 1.57079632, 'diffractionAngle'],

    # After_HRM1_HRM2: drift
    ['op_After_HRM1_HRM2_L', 'f', 0.07999999999998408, 'length'],

    # HRM2: crystal
    ['op_HRM2_hfn', 's', '', 'heightProfileFile'],
    ['op_HRM2_dim', 's', 'x', 'orientation'],
    ['op_HRM2_d_sp', 'f', 0.9600687344008111, 'dSpacing'],
    ['op_HRM2_psi0r', 'f', -1.0873035035585694e-05, 'psi0r'],
    ['op_HRM2_psi0i', 'f', 1.8438837339536554e-07, 'psi0i'],
    ['op_HRM2_psiHr', 'f', -4.181686438547451e-06, 'psiHr'],
    ['op_HRM2_psiHi', 'f', 1.6100412693351052e-07, 'psiHi'],
    ['op_HRM2_psiHBr', 'f', -4.181686438547451e-06, 'psiHBr'],
    ['op_HRM2_psiHBi', 'f', 1.6100412693351052e-07, 'psiHBi'],
    ['op_HRM2_tc', 'f', 0.01, 'crystalThickness'],
    ['op_HRM2_uc', 'f', 1, 'useCase'],
    ['op_HRM2_ang_as', 'f', -0.51487212, 'asymmetryAngle'],
    ['op_HRM2_nvx', 'f', 0.30193330638751636, 'nvx'],
    ['op_HRM2_nvy', 'f', 2.051605629178758e-09, 'nvy'],
    ['op_HRM2_nvz', 'f', -0.9533290504825195, 'nvz'],
    ['op_HRM2_tvx', 'f', 0.9533290504825195, 'tvx'],
    ['op_HRM2_tvy', 'f', 6.477772425408857e-09, 'tvy'],
    ['op_HRM2_ang', 'f', 1.2640763693590809, 'grazingAngle'],
    ['op_HRM2_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_HRM2_energy', 'f', 9481.0, 'energy'],
    ['op_HRM2_diffractionAngle', 'f', -1.57079632, 'diffractionAngle'],

    # After_HRM2_HRM_F1: drift
    ['op_After_HRM2_HRM_F1_L', 'f', 10.0, 'length'],

    # HRM_F1: crl
    ['op_HRM_F1_foc_plane', 'f', 1, 'focalPlane'],
    ['op_HRM_F1_delta', 'f', 3.791135e-06, 'refractiveIndex'],
    ['op_HRM_F1_atten_len', 'f', 0.008387, 'attenuationLength'],
    ['op_HRM_F1_shape', 'f', 1, 'shape'],
    ['op_HRM_F1_apert_h', 'f', 0.0015, 'horizontalApertureSize'],
    ['op_HRM_F1_apert_v', 'f', 0.0015, 'verticalApertureSize'],
    ['op_HRM_F1_r_min', 'f', 0.0002, 'tipRadius'],
    ['op_HRM_F1_wall_thick', 'f', 5e-05, 'tipWallThickness'],
    ['op_HRM_F1_x', 'f', 0.0, 'horizontalOffset'],
    ['op_HRM_F1_y', 'f', 0.0, 'verticalOffset'],
    ['op_HRM_F1_n', 'i', 1, 'numberOfLenses'],

    # After_HRM_F1_Before_WDS: drift
    ['op_After_HRM_F1_Before_WDS_L', 'f', 10.0, 'length'],

    # WDS: aperture
    ['op_WDS_shape', 's', 'r', 'shape'],
    ['op_WDS_Dx', 'f', 0.0002, 'horizontalSize'],
    ['op_WDS_Dy', 'f', 0.01, 'verticalSize'],
    ['op_WDS_x', 'f', 0.0, 'horizontalOffset'],
    ['op_WDS_y', 'f', 0.0, 'verticalOffset'],

    # After_WDS_HRM_F2: drift
    ['op_After_WDS_HRM_F2_L', 'f', 10.0, 'length'],

    # HRM_F2: crl
    ['op_HRM_F2_foc_plane', 'f', 1, 'focalPlane'],
    ['op_HRM_F2_delta', 'f', 3.791135e-06, 'refractiveIndex'],
    ['op_HRM_F2_atten_len', 'f', 0.008387, 'attenuationLength'],
    ['op_HRM_F2_shape', 'f', 1, 'shape'],
    ['op_HRM_F2_apert_h', 'f', 0.0015, 'horizontalApertureSize'],
    ['op_HRM_F2_apert_v', 'f', 0.0015, 'verticalApertureSize'],
    ['op_HRM_F2_r_min', 'f', 0.0002, 'tipRadius'],
    ['op_HRM_F2_wall_thick', 'f', 5e-05, 'tipWallThickness'],
    ['op_HRM_F2_x', 'f', 0.0, 'horizontalOffset'],
    ['op_HRM_F2_y', 'f', 0.0, 'verticalOffset'],
    ['op_HRM_F2_n', 'i', 1, 'numberOfLenses'],

    # After_HRM_F2_Before_HRM3: drift
    ['op_After_HRM_F2_Before_HRM3_L', 'f', 10.0, 'length'],

    # HRM3: crystal
    ['op_HRM3_hfn', 's', '', 'heightProfileFile'],
    ['op_HRM3_dim', 's', 'x', 'orientation'],
    ['op_HRM3_d_sp', 'f', 0.9600687344008111, 'dSpacing'],
    ['op_HRM3_psi0r', 'f', -1.0873035035585694e-05, 'psi0r'],
    ['op_HRM3_psi0i', 'f', 1.8438837339536554e-07, 'psi0i'],
    ['op_HRM3_psiHr', 'f', -4.181686438547451e-06, 'psiHr'],
    ['op_HRM3_psiHi', 'f', 1.6100412693351052e-07, 'psiHi'],
    ['op_HRM3_psiHBr', 'f', -4.181686438547451e-06, 'psiHBr'],
    ['op_HRM3_psiHBi', 'f', 1.6100412693351052e-07, 'psiHBi'],
    ['op_HRM3_tc', 'f', 0.01, 'crystalThickness'],
    ['op_HRM3_uc', 'f', 1, 'useCase'],
    ['op_HRM3_ang_as', 'f', 0.51487212, 'asymmetryAngle'],
    ['op_HRM3_nvx', 'f', 0.972664744212603, 'nvx'],
    ['op_HRM3_nvy', 'f', 6.609156467054803e-09, 'nvy'],
    ['op_HRM3_nvz', 'f', -0.23221390002717662, 'nvz'],
    ['op_HRM3_tvx', 'f', 0.23221390002717662, 'tvx'],
    ['op_HRM3_tvy', 'f', 1.5778694645163085e-09, 'tvy'],
    ['op_HRM3_ang', 'f', 0.23435318504297908, 'grazingAngle'],
    ['op_HRM3_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_HRM3_energy', 'f', 9481.0, 'energy'],
    ['op_HRM3_diffractionAngle', 'f', -1.57079632, 'diffractionAngle'],

    # After_HRM3_HRM4: drift
    ['op_After_HRM3_HRM4_L', 'f', 0.07999999999998408, 'length'],

    # HRM4: crystal
    ['op_HRM4_hfn', 's', '', 'heightProfileFile'],
    ['op_HRM4_dim', 's', 'x', 'orientation'],
    ['op_HRM4_d_sp', 'f', 0.9600687344008111, 'dSpacing'],
    ['op_HRM4_psi0r', 'f', -1.0873035035585694e-05, 'psi0r'],
    ['op_HRM4_psi0i', 'f', 1.8438837339536554e-07, 'psi0i'],
    ['op_HRM4_psiHr', 'f', -4.181686438547451e-06, 'psiHr'],
    ['op_HRM4_psiHi', 'f', 1.6100412693351052e-07, 'psiHi'],
    ['op_HRM4_psiHBr', 'f', -4.181686438547451e-06, 'psiHBr'],
    ['op_HRM4_psiHBi', 'f', 1.6100412693351052e-07, 'psiHBi'],
    ['op_HRM4_tc', 'f', 0.01, 'crystalThickness'],
    ['op_HRM4_uc', 'f', 1, 'useCase'],
    ['op_HRM4_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_HRM4_nvx', 'f', -0.7322282430733594, 'nvx'],
    ['op_HRM4_nvy', 'f', 4.975415277322606e-09, 'nvy'],
    ['op_HRM4_nvz', 'f', -0.6810593219725439, 'nvz'],
    ['op_HRM4_tvx', 'f', -0.6810593219725439, 'tvx'],
    ['op_HRM4_tvy', 'f', 4.627727743855522e-09, 'tvy'],
    ['op_HRM4_ang', 'f', 0.7492083731847909, 'grazingAngle'],
    ['op_HRM4_amp_coef', 'f', 1.0, 'heightAmplification'],
    ['op_HRM4_energy', 'f', 9481.0, 'energy'],
    ['op_HRM4_diffractionAngle', 'f', 1.57079632, 'diffractionAngle'],

    # After_HRM4_After_HRM: drift
    ['op_After_HRM4_After_HRM_L', 'f', 10.0, 'length'],

#---Propagation parameters
    ['op_MR1L3_pp', 'f',                    [0, 0, 1.0, 0, 0, 1.0, 4.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'MR1L3'],
    ['op_MR1L3_MR2L3_pp', 'f',              [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'MR1L3_MR2L3'],
    ['op_MR2L3_pp', 'f',                    [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'MR2L3'],
    ['op_MR2L3_Before_HHLM_pp', 'f',        [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'MR2L3_Before_HHLM'],
    ['op_C1_pp', 'f',                       [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, -0.640439232226, 4.352e-09, 0.768008847492, 0.768008847492, 1.576e-09], 'C1'],
    ['op_After_C1_C2_pp', 'f',              [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_C1_C2'],
    ['op_C2_pp', 'f',                       [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.997381741513, 6.777e-09, 0.072316399909, 0.072316399909, -6.304e-09], 'C2'],
    ['op_After_C2_C3_pp', 'f',              [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_C2_C3'],
    ['op_C3_pp', 'f',                       [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, -0.997381741513, 6.777e-09, 0.072316399909, 0.072316399909, 6.304e-09], 'C3'],
    ['op_After_C3_C4_pp', 'f',              [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_C3_C4'],
    ['op_C4_pp', 'f',                       [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.640439135646, 4.352e-09, 0.768008928029, 0.768008928029, -1.576e-09], 'C4'],
    ['op_After_C4_HRM1_pp', 'f',            [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_C4_HRM1'],
    ['op_HRM1_pp', 'f',                     [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, -0.997381741513, 6.777e-09, 0.072316399909, 0.072316399909, 6.304e-09], 'HRM1'],
    ['op_After_HRM1_HRM2_pp', 'f',          [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_HRM1_HRM2'],
    ['op_HRM2_pp', 'f',                     [0, 0, 1.0, 0, 0, 3.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.997382667547, 6.777e-09, 0.072303626994, 0.072303626994, -6.304e-09], 'HRM2'],
    ['op_After_HRM2_HRM_F1_pp', 'f',        [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_HRM2_HRM_F1'],
    ['op_HRM_F1_pp', 'f',                   [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HRM_F1'],
    ['op_After_HRM_F1_Before_WDS_pp', 'f',  [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_HRM_F1_Before_WDS'],
    ['op_WDS_pp', 'f',                      [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'WDS'],
    ['op_After_WDS_HRM_F2_pp', 'f',         [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_WDS_HRM_F2'],
    ['op_HRM_F2_pp', 'f',                   [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HRM_F2'],
    ['op_After_HRM_F2_Before_HRM3_pp', 'f', [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_HRM_F2_Before_HRM3'],
    ['op_HRM3_pp', 'f',                     [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.99738266769, 6.777e-09, 0.072303625018, 0.072303625018, -6.304e-09], 'HRM3'],
    ['op_After_HRM3_HRM4_pp', 'f',          [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_HRM3_HRM4'],
    ['op_HRM4_pp', 'f',                     [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, -0.997381741513, 6.777e-09, 0.072316399909, 0.072316399909, 6.304e-09], 'HRM4'],
    ['op_After_HRM4_After_HRM_pp', 'f',     [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_HRM4_After_HRM'],
    ['op_fin_pp', 'f',                      [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'final post-propagation (resize) parameters'],

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


def main():
    v = srwl_bl.srwl_uti_parse_options(varParam, use_sys_argv=True)
    op = set_optics(v)
    v.si = True
    v.si_pl = 'xy'
    v.ws = True
    v.ws_pl = 'xy'
    mag = None
    if v.rs_type == 'm':
        mag = srwlib.SRWLMagFldC()
        mag.arXc.append(0)
        mag.arYc.append(0)
        mag.arMagFld.append(srwlib.SRWLMagFldM(v.mp_field, v.mp_order, v.mp_distribution, v.mp_len))
        mag.arZc.append(v.mp_zc)
    srwl_bl.SRWLBeamline(_name=v.name, _mag_approx=mag).calc_all(v, op)

main()
