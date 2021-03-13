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

from time import *
from copy import *
from array import *
from uti_plot import *

#------------------------------------------------------------------------------
def set_optics_before_split(v=None):
    el = []
    pp = []
    names = ['Aperture', 'Aperture_Split']
    for el_name in names:
        if el_name == 'Aperture':
            # Aperture: aperture 399.3m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_Aperture_shape,
                _ap_or_ob='a',
                _Dx=v.op_Aperture_Dx,
                _Dy=v.op_Aperture_Dy,
                _x=v.op_Aperture_x,
                _y=v.op_Aperture_y,
            ))
            pp.append(v.op_Aperture_pp)
        elif el_name == 'Aperture_Split':
            # Aperture_Split: drift 399.3m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Aperture_Split_L,
                #_L=v.op_Aperture_Split_CC_L,
            ))
            pp.append(v.op_Aperture_Split_pp)
            #pp.append(v.op_Aperture_Split_CC_pp)
    return srwlib.SRWLOptC(el, pp)
             
#------------------------------------------------------------------------------
def set_optics_cc(v=None):
    el = []
    pp = []
    names = ['Split_CC', 'CC1_1', 'CC1_1_CC1_2', 'CC1_2', 'CC1_2_CC2_1', 'CC2_1', 'CC2_1_CC2_2', 'CC2_2', 'Recomb_CC']
    #names = ['Split_CC', 'CC1_1', 'CC1_1_CC1_2', 'CC1_2', 'Watchpoint', 'Watchpoint_CC2_1', 'CC2_1', 'CC2_1_CC2_2', 'CC2_2', 'Recomb_CC', 'Recomb_CC_Watchpoint2', 'Watchpoint2', 'CRL', 'Watchpoint3', 'Watchpoint3_Monitor', 'Monitor']
    for el_name in names:
        if el_name == 'Split_CC':
            # Split_CC: aperture 400.1m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_Split_CC_shape,
                _ap_or_ob='a',
                _Dx=v.op_Split_CC_Dx,
                _Dy=v.op_Split_CC_Dy,
                _x=v.op_Split_CC_x,
                _y=v.op_Split_CC_y,
            ))
            pp.append(v.op_Split_CC_pp)
        elif el_name == 'CC1_1':
            # CC1_1: crystal 400.1m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_CC1_1_d_sp,
                _psi0r=v.op_CC1_1_psi0r,
                _psi0i=v.op_CC1_1_psi0i,
                _psi_hr=v.op_CC1_1_psiHr,
                _psi_hi=v.op_CC1_1_psiHi,
                _psi_hbr=v.op_CC1_1_psiHBr,
                _psi_hbi=v.op_CC1_1_psiHBi,
                _tc=v.op_CC1_1_tc,
                _ang_as=v.op_CC1_1_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_CC1_1_nvx,
                _nvy=v.op_CC1_1_nvy,
                _nvz=v.op_CC1_1_nvz,
                _tvx=v.op_CC1_1_tvx,
                _tvy=v.op_CC1_1_tvy,
            )
            el.append(crystal)
            pp.append(v.op_CC1_1_pp)
            
        elif el_name == 'CC1_1_CC1_2':
            # CC1_1_CC1_2: drift 400.1m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CC1_1_CC1_2_L,
            ))
            pp.append(v.op_CC1_1_CC1_2_pp)
        elif el_name == 'CC1_2':
            # CC1_2: crystal 400.1127m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_CC1_2_d_sp,
                _psi0r=v.op_CC1_2_psi0r,
                _psi0i=v.op_CC1_2_psi0i,
                _psi_hr=v.op_CC1_2_psiHr,
                _psi_hi=v.op_CC1_2_psiHi,
                _psi_hbr=v.op_CC1_2_psiHBr,
                _psi_hbi=v.op_CC1_2_psiHBi,
                _tc=v.op_CC1_2_tc,
                _ang_as=v.op_CC1_2_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_CC1_2_nvx,
                _nvy=v.op_CC1_2_nvy,
                _nvz=v.op_CC1_2_nvz,
                _tvx=v.op_CC1_2_tvx,
                _tvy=v.op_CC1_2_tvy,
            )
            el.append(crystal)
            pp.append(v.op_CC1_2_pp)
            
        elif el_name == 'CC1_2_CC2_1':
        #elif el_name == 'Watchpoint_CC2_1':
            # Watchpoint_CC2_1: drift 400.1127m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CC1_2_CC2_1_L,
            ))
            pp.append(v.op_CC1_2_CC2_1_pp)
            
            if(v.op_CC_ang_y != 0.):
                el.append(srwlib.SRWLOptAng(_ang_y=v.op_CC_ang_y))
                pp.append(v.op_CC_ang_pp)
            
        elif el_name == 'CC2_1':
            # CC2_1: crystal 400.7731m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_CC2_1_d_sp,
                _psi0r=v.op_CC2_1_psi0r,
                _psi0i=v.op_CC2_1_psi0i,
                _psi_hr=v.op_CC2_1_psiHr,
                _psi_hi=v.op_CC2_1_psiHi,
                _psi_hbr=v.op_CC2_1_psiHBr,
                _psi_hbi=v.op_CC2_1_psiHBi,
                _tc=v.op_CC2_1_tc,
                _ang_as=v.op_CC2_1_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_CC2_1_nvx,
                _nvy=v.op_CC2_1_nvy,
                _nvz=v.op_CC2_1_nvz,
                _tvx=v.op_CC2_1_tvx,
                _tvy=v.op_CC2_1_tvy,
            )
            el.append(crystal)
            pp.append(v.op_CC2_1_pp)
            
        elif el_name == 'CC2_1_CC2_2':
            # CC2_1_CC2_2: drift 400.7731m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CC2_1_CC2_2_L,
            ))
            pp.append(v.op_CC2_1_CC2_2_pp)
        elif el_name == 'CC2_2':
            # CC2_2: crystal 400.7858m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_CC2_2_d_sp,
                _psi0r=v.op_CC2_2_psi0r,
                _psi0i=v.op_CC2_2_psi0i,
                _psi_hr=v.op_CC2_2_psiHr,
                _psi_hi=v.op_CC2_2_psiHi,
                _psi_hbr=v.op_CC2_2_psiHBr,
                _psi_hbi=v.op_CC2_2_psiHBi,
                _tc=v.op_CC2_2_tc,
                _ang_as=v.op_CC2_2_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_CC2_2_nvx,
                _nvy=v.op_CC2_2_nvy,
                _nvz=v.op_CC2_2_nvz,
                _tvx=v.op_CC2_2_tvx,
                _tvy=v.op_CC2_2_tvy,
            )
            el.append(crystal)
            pp.append(v.op_CC2_2_pp)
            
        elif el_name == 'Recomb_CC':
            # Recomb_CC: aperture 400.7858m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_Recomb_CC_shape,
                _ap_or_ob='a',
                _Dx=v.op_Recomb_CC_Dx,
                _Dy=v.op_Recomb_CC_Dy,
                _x=v.op_Recomb_CC_x,
                _y=v.op_Recomb_CC_y,
            ))
            pp.append(v.op_Recomb_CC_pp)
            
    return srwlib.SRWLOptC(el, pp)

#------------------------------------------------------------------------------
def set_optics_vcc(v=None):
    el = []
    pp = []
    names = ['Split_VCC', 'Split_VCC_VCC1_1', 'VCC1_1', 'VCC1_1_VCC1_2', 'VCC1_2', 'VCC1_2_VCC2_1', 'VCC2_1', 'VCC2_1_VCC2_2', 'VCC2_2', 'VCC2_2_Recomb_VCC', 'Recomb_VCC']
    #names = ['Split_VCC', 'Split_VCC_VCC1_1', 'VCC1_1', 'VCC1_1_VCC1_2', 'VCC1_2', 'Watchpoint', 'Watchpoint_VCC2_1', 'VCC2_1', 'VCC2_1_VCC2_2', 'VCC2_2', 'VCC2_2_Recomb_VCC', 'Recomb_VCC', 'Recomb_VCC_Watchpoint2', 'Watchpoint2', 'CRL', 'Watchpoint3', 'Watchpoint3_Monitor', 'Monitor']
    for el_name in names:
        if el_name == 'Split_VCC':
            # Split_VCC: aperture 400.1m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_Split_VCC_shape,
                _ap_or_ob='a',
                _Dx=v.op_Split_VCC_Dx,
                _Dy=v.op_Split_VCC_Dy,
                _x=v.op_Split_VCC_x,
                _y=v.op_Split_VCC_y,
            ))
            pp.append(v.op_Split_VCC_pp)
        elif el_name == 'Split_VCC_VCC1_1':
            # Split_VCC_VCC1_1: drift 400.1m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Split_VCC_VCC1_1_L,
            ))
            pp.append(v.op_Split_VCC_VCC1_1_pp)
        elif el_name == 'VCC1_1':
            # VCC1_1: crystal 400.3m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_VCC1_1_d_sp,
                _psi0r=v.op_VCC1_1_psi0r,
                _psi0i=v.op_VCC1_1_psi0i,
                _psi_hr=v.op_VCC1_1_psiHr,
                _psi_hi=v.op_VCC1_1_psiHi,
                _psi_hbr=v.op_VCC1_1_psiHBr,
                _psi_hbi=v.op_VCC1_1_psiHBi,
                _tc=v.op_VCC1_1_tc,
                _ang_as=v.op_VCC1_1_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_VCC1_1_nvx,
                _nvy=v.op_VCC1_1_nvy,
                _nvz=v.op_VCC1_1_nvz,
                _tvx=v.op_VCC1_1_tvx,
                _tvy=v.op_VCC1_1_tvy,
            )
            el.append(crystal)
            pp.append(v.op_VCC1_1_pp)
            
        elif el_name == 'VCC1_1_VCC1_2':
            # VCC1_1_VCC1_2: drift 400.3m
            el.append(srwlib.SRWLOptD(
                _L=v.op_VCC1_1_VCC1_2_L,
            ))
            pp.append(v.op_VCC1_1_VCC1_2_pp)
        elif el_name == 'VCC1_2':
            # VCC1_2: crystal 400.3127m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_VCC1_2_d_sp,
                _psi0r=v.op_VCC1_2_psi0r,
                _psi0i=v.op_VCC1_2_psi0i,
                _psi_hr=v.op_VCC1_2_psiHr,
                _psi_hi=v.op_VCC1_2_psiHi,
                _psi_hbr=v.op_VCC1_2_psiHBr,
                _psi_hbi=v.op_VCC1_2_psiHBi,
                _tc=v.op_VCC1_2_tc,
                _ang_as=v.op_VCC1_2_ang_as,
            )
            
            pi_rad = 3.14159265358979
            orientDataVCC1_2 = crystal.find_orient(v.gbm_ave, _ang_dif_pl=pi_rad/2)
            orientVCC1_2 = orientDataVCC1_2[0] #VCC1_2 crystal orientation found
            tVCC1_2 = orientVCC1_2[0]; nVCC1_2 = orientVCC1_2[2] #Tangential and Normal vectors to crystal surface
            #print('VCC1_2 crystal orientation found:')
            #print('t =', tVCC1_2, 's =', orientVCC1_2[1], 'n =', nVCC1_2)
            #print('VCC1_2 crystal orientation defined:')
            #print('[tx,ty] =', [v.op_VCC1_2_tvx, v.op_VCC1_2_tvy], 'n =', [v.op_VCC1_2_nvx, v.op_VCC1_2_nvy, v.op_VCC1_2_nvz])            

            #Set orientation of VCC1_2 crystal:
            crystal.set_orient(
                _nvx=nVCC1_2[0], 
                _nvy=nVCC1_2[1], 
                _nvz=nVCC1_2[2], 
                _tvx=tVCC1_2[0], 
                _tvy=tVCC1_2[1],
            )
            #crystal.set_orient(
            #    _nvx=v.op_VCC1_2_nvx,
            #    _nvy=v.op_VCC1_2_nvy,
            #    _nvz=v.op_VCC1_2_nvz,
            #    _tvx=v.op_VCC1_2_tvx,
            #    _tvy=v.op_VCC1_2_tvy,
            #)
            el.append(crystal)
            
            orientOutFrVCC1_2 = orientDataVCC1_2[1] #Orientation of the Outgoing beam frame being found
            exVCC1_2 = orientOutFrVCC1_2[0]; eyVCC1_2 = orientOutFrVCC1_2[1]; ezVCC1_2 = orientOutFrVCC1_2[2] #Horizontal, Vertical and Longitudinal base vectors of the Output beam frame
            #print('VCC1_2 crystal Outgoing beam frame:'); 
            #print('ex =', exVCC1_2, 'ey =', eyVCC1_2, 'ez =', ezVCC1_2)
            v.op_VCC1_2_pp[12] = ezVCC1_2[0]; v.op_VCC1_2_pp[13] = ezVCC1_2[1]; v.op_VCC1_2_pp[14] = ezVCC1_2[2]
            v.op_VCC1_2_pp[15] = exVCC1_2[0]; v.op_VCC1_2_pp[16] = exVCC1_2[1]
            #[12]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Horizontal Coordinate
            #[13]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Vertical Coordinate
            #[14]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Longitudinal Coordinate
            #[15]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Horizontal Coordinate
            #[16]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Vertical Coordinate
            pp.append(v.op_VCC1_2_pp)
            
        elif el_name == 'VCC1_2_VCC2_1':
        #elif el_name == 'Watchpoint_VCC2_1':
            # Watchpoint_VCC2_1: drift 400.3127m
            el.append(srwlib.SRWLOptD(
                _L=v.op_VCC1_2_VCC2_1_L,
                #_L=v.op_Watchpoint_VCC2_1_L,
            ))
            pp.append(v.op_VCC1_2_VCC2_1_pp)

            if(v.op_VCC1_2_VCC2_1_dt != 0.):
                el.append(srwlib.SRWLOptD(_L=(v.op_VCC1_2_VCC2_1_dt*2.99792458e+08), _treat=1))
                pp.append(v.op_VCC1_2_VCC2_1_dt_pp)

            if((v.op_VCC_shift_x != 0.) or (v.op_VCC_shift_y != 0.)):
                el.append(srwlib.SRWLOptShift(_shift_x=v.op_VCC_shift_x, _shift_y=v.op_VCC_shift_y))
                pp.append(v.op_VCC_shift_pp)
                
            #pp.append(v.op_Watchpoint_VCC2_1_pp)
        elif el_name == 'VCC2_1':
            # VCC2_1: crystal 400.5731m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_VCC2_1_d_sp,
                _psi0r=v.op_VCC2_1_psi0r,
                _psi0i=v.op_VCC2_1_psi0i,
                _psi_hr=v.op_VCC2_1_psiHr,
                _psi_hi=v.op_VCC2_1_psiHi,
                _psi_hbr=v.op_VCC2_1_psiHBr,
                _psi_hbi=v.op_VCC2_1_psiHBi,
                _tc=v.op_VCC2_1_tc,
                _ang_as=v.op_VCC2_1_ang_as,
            )
            crystal.set_orient(
                _nvx=v.op_VCC2_1_nvx,
                _nvy=v.op_VCC2_1_nvy,
                _nvz=v.op_VCC2_1_nvz,
                _tvx=v.op_VCC2_1_tvx,
                _tvy=v.op_VCC2_1_tvy,
            )
            el.append(crystal)
            pp.append(v.op_VCC2_1_pp)
            
        elif el_name == 'VCC2_1_VCC2_2':
            # VCC2_1_VCC2_2: drift 400.5731m
            el.append(srwlib.SRWLOptD(
                _L=v.op_VCC2_1_VCC2_2_L,
            ))
            pp.append(v.op_VCC2_1_VCC2_2_pp)
        elif el_name == 'VCC2_2':
            # VCC2_2: crystal 400.5858m
            crystal = srwlib.SRWLOptCryst(
                _d_sp=v.op_VCC2_2_d_sp,
                _psi0r=v.op_VCC2_2_psi0r,
                _psi0i=v.op_VCC2_2_psi0i,
                _psi_hr=v.op_VCC2_2_psiHr,
                _psi_hi=v.op_VCC2_2_psiHi,
                _psi_hbr=v.op_VCC2_2_psiHBr,
                _psi_hbi=v.op_VCC2_2_psiHBi,
                _tc=v.op_VCC2_2_tc,
                _ang_as=v.op_VCC2_2_ang_as,
            )
            
            orientDataVCC2_2 = crystal.find_orient(v.gbm_ave, _ang_dif_pl=-pi_rad/2)
            orientVCC2_2 = orientDataVCC2_2[0] #VCC2_2 crystal orientation found
            tVCC2_2 = orientVCC2_2[0]; nVCC2_2 = orientVCC2_2[2] #Tangential and Normal vectors to crystal surface
            #print('VCC2_2 crystal orientation found:')
            #print('t =', tVCC2_2, 's =', orientVCC2_2[1], 'n =', nVCC2_2)
            #print('VCC2_2 crystal orientation defined:')
            #print('[tx,ty] =', [v.op_VCC2_2_tvx, v.op_VCC2_2_tvy], 'n =', [v.op_VCC2_2_nvx, v.op_VCC2_2_nvy, v.op_VCC2_2_nvz])            
            #Set orientation of VCC2_2 crystal:
            crystal.set_orient(
                _nvx=nVCC2_2[0], 
                _nvy=nVCC2_2[1], 
                _nvz=nVCC2_2[2], 
                _tvx=tVCC2_2[0], 
                _tvy=tVCC2_2[1],
            )
            #crystal.set_orient(
            #    _nvx=v.op_VCC2_2_nvx,
            #    _nvy=v.op_VCC2_2_nvy,
            #    _nvz=v.op_VCC2_2_nvz,
            #    _tvx=v.op_VCC2_2_tvx,
            #    _tvy=v.op_VCC2_2_tvy,
            #)
            el.append(crystal)
            
            orientOutFrVCC2_2 = orientDataVCC2_2[1] #Orientation of the Outgoing beam frame being found
            exVCC2_2 = orientOutFrVCC2_2[0]; eyVCC2_2 = orientOutFrVCC2_2[1]; ezVCC2_2 = orientOutFrVCC2_2[2] #Horizontal, Vertical and Longitudinal base vectors of the Output beam frame
            #print('VCC2_2 crystal Outgoing beam frame:'); 
            #print('ex =', exVCC2_2, 'ey =', eyVCC2_2, 'ez =', ezVCC2_2)
            v.op_VCC2_2_pp[12] = ezVCC2_2[0]; v.op_VCC2_2_pp[13] = ezVCC2_2[1]; v.op_VCC2_2_pp[14] = ezVCC2_2[2]
            v.op_VCC2_2_pp[15] = exVCC2_2[0]; v.op_VCC2_2_pp[16] = exVCC2_2[1]
            #[12]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Horizontal Coordinate
            #[13]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Vertical Coordinate
            #[14]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Longitudinal Coordinate
            #[15]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Horizontal Coordinate
            #[16]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Vertical Coordinate
            pp.append(v.op_VCC2_2_pp)
            
        elif el_name == 'VCC2_2_Recomb_VCC':
            # VCC2_2_Recomb_VCC: drift 400.5858m
            el.append(srwlib.SRWLOptD(
                _L=v.op_VCC2_2_Recomb_VCC_L,
            ))
            pp.append(v.op_VCC2_2_Recomb_VCC_pp)
        elif el_name == 'Recomb_VCC':
            # Recomb_VCC: aperture 400.7858m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_Recomb_VCC_shape,
                _ap_or_ob='a',
                _Dx=v.op_Recomb_VCC_Dx,
                _Dy=v.op_Recomb_VCC_Dy,
                _x=v.op_Recomb_VCC_x,
                _y=v.op_Recomb_VCC_y,
            ))
            pp.append(v.op_Recomb_VCC_pp)
            
    return srwlib.SRWLOptC(el, pp)

#------------------------------------------------------------------------------
def set_optics_after_recomb(v=None, _use_crl=True):
    el = []
    pp = []
    names = ['Recomb_CRL', 'CRL1', 'CRL2', 'CRL3', 'CRL_Monitor']
    if(not _use_crl):
        names = ['Recomb_CRL']
        v.op_Recomb_CRL_L += v.op_CRL_Monitor_L
    
    for el_name in names:
            
        if el_name == 'Recomb_CRL':
        #if el_name == 'Recomb_VCC_Watchpoint2':
            # Recomb_VCC_Watchpoint2: drift 400.7858m
            el.append(srwlib.SRWLOptD(
                _L=v.op_Recomb_CRL_L,
            ))
            pp.append(v.op_Recomb_CRL_pp)
            #pp.append(v.op_Recomb_VCC_Watchpoint2_pp)

        elif el_name == 'CRL1':
            # CRL1: crl 406.3858m
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
        elif el_name == 'CRL2':
            # CRL2: crl 406.3858m
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
        elif el_name == 'CRL3':
            # CRL3: crl 406.3858m
            el.append(srwlib.srwl_opt_setup_CRL(
                _foc_plane=v.op_CRL3_foc_plane,
                _delta=v.op_CRL3_delta,
                _atten_len=v.op_CRL3_atten_len,
                _shape=v.op_CRL3_shape,
                _apert_h=v.op_CRL3_apert_h,
                _apert_v=v.op_CRL3_apert_v,
                _r_min=v.op_CRL3_r_min,
                _n=v.op_CRL3_n,
                _wall_thick=v.op_CRL3_wall_thick,
                _xc=v.op_CRL3_x,
                _yc=v.op_CRL3_y,
            ))
            pp.append(v.op_CRL3_pp)
        elif el_name == 'CRL_Monitor':
            # CRL_Monitor: drift 406.3858m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CRL_Monitor_L,
            ))
            pp.append(v.op_CRL_Monitor_pp)
    pp.append(v.op_fin_pp)
    return srwlib.SRWLOptC(el, pp)

#------------------------------------------------------------------------------
varParam = srwl_bl.srwl_uti_ext_options([
    ['name', 's', 'Split-Delay beased on Channel-Cut Crystals CC Branch 01', 'simulation name'],

#---Data Folder
    ['fdir', 's', '', 'folder (directory) name for reading-in input and saving output data files'],

    ['gbm_x', 'f', 0.0, 'average horizontal coordinates of waist [m]'],
    ['gbm_y', 'f', 0.0, 'average vertical coordinates of waist [m]'],
    ['gbm_z', 'f', 0.0, 'average longitudinal coordinate of waist [m]'],
    ['gbm_xp', 'f', 0.0, 'average horizontal angle at waist [rad]'],
    ['gbm_yp', 'f', 0.0, 'average verical angle at waist [rad]'],
    ['gbm_ave', 'f', 9500.0, 'average photon energy [eV]'],
    ['gbm_pen', 'f', 0.001, 'energy per pulse [J]'],
    ['gbm_rep', 'f', 1, 'rep. rate [Hz]'],
    ['gbm_pol', 'f', 1, 'polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left'],
    ['gbm_sx', 'f', 8.221943e-06, 'rms beam size vs horizontal position [m] at waist (for intensity)'],
    ['gbm_sy', 'f', 8.221943e-06, 'rms beam size vs vertical position [m] at waist (for intensity)'],
    ['gbm_st', 'f', 1e-14, 'rms pulse duration [s] (for intensity)'],
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

    #['w_e', 'f', 9500.0, 'photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    #['w_ef', 'f', -1.0, 'final photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    #['w_ne', 'i', 1, 'number of points vs photon energy for calculation of intensity distribution'],
    ['w_e', 'f', -4e-14, 'photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ef', 'f', 4e-14, 'final photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ne', 'i', 100, 'number of points vs photon energy or time for calculation of intensity distribution'],

    ['w_x', 'f', 0.0, 'central horizontal position [m] for calculation of intensity distribution'],
    ['w_rx', 'f', 0.004, 'range of horizontal position [m] for calculation of intensity distribution'],
    ['w_nx', 'i', 100, 'number of points vs horizontal position for calculation of intensity distribution'],
    ['w_y', 'f', 0.0, 'central vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ry', 'f', 0.004, 'range of vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ny', 'i', 100, 'number of points vs vertical position for calculation of intensity distribution'],
    #['w_smpf', 'f', 0.3, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_smpf', 'f', -1, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_meth', 'i', 2, 'method to use for calculation of intensity distribution vs horizontal and vertical position: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['w_prec', 'f', 0.01, 'relative precision for calculation of intensity distribution vs horizontal and vertical position'],

    ['w_ft', 's', 't', 'presentation/domain: "f"- frequency (photon energy), "t"- time'],

    ['w_u', 'i', 2, 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['si_pol', 'i', 6, 'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['si_type', 'i', 0, 'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],
    #['w_mag', 'i', 1, 'magnetic field to be used for calculation of intensity distribution vs horizontal and vertical position: 1- approximate, 2- accurate'],

    ['si_fn', 's', 'res_int_se.dat', 'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],
    ['si_pl', 's', '', 'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
    ['ws_fni', 's', 'res_int_pr_se.dat', 'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],
    ['ws_pl', 's', '', 'plot the resulting intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    #to add options
    ['op_r', 'f', 399.3, 'longitudinal position of the first optical element [m]'],
    #['op_r', 'f', 400.1, 'longitudinal position of the first optical element [m]'],
    
    #['op_fno', 's', 'spit_delay_bl_orientation.dat', 'file name for saving orientations of optical elements in the lab frame'],


    # Former appParam:
    ['rs_type', 's', 'g', 'source type, (u) idealized undulator, (t), tabulated undulator, (m) multipole, (g) gaussian beam'],

#---Beamline optics:

    ##Before Split
    # Aperture: aperture
    ['op_Aperture_shape', 's', 'r', 'shape'],
    ['op_Aperture_Dx', 'f', 0.001, 'horizontalSize'],
    ['op_Aperture_Dy', 'f', 0.001, 'verticalSize'],
    ['op_Aperture_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Aperture_y', 'f', 0.0, 'verticalOffset'],

    # Aperture_Split: drift
    ['op_Aperture_Split_L', 'f', 0.8, 'length'],

    ##CC Branch
    # Split_CC: aperture
    ['op_Split_CC_shape', 's', 'r', 'shape'],
    ['op_Split_CC_Dx', 'f', 0.01, 'horizontalSize'],
    ['op_Split_CC_Dy', 'f', 0.004, 'verticalSize'],
    ['op_Split_CC_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Split_CC_y', 'f', 0.002, 'verticalOffset'],

    # CC1_1: crystal
    ['op_CC1_1_hfn', 's', '', 'heightProfileFile'],
    ['op_CC1_1_dim', 's', 'x', 'orientation'],
    ['op_CC1_1_d_sp', 'f', 1.9201374688, 'dSpacing'],
    ['op_CC1_1_psi0r', 'f', -1.08291732223e-05, 'psi0r'],
    ['op_CC1_1_psi0i', 'f', 1.82929605452e-07, 'psi0i'],
    ['op_CC1_1_psiHr', 'f', -6.58335380584e-06, 'psiHr'],
    ['op_CC1_1_psiHi', 'f', 1.76831595902e-07, 'psiHi'],
    ['op_CC1_1_psiHBr', 'f', -6.58335380584e-06, 'psiHBr'],
    ['op_CC1_1_psiHBi', 'f', 1.76831595902e-07, 'psiHBi'],
    ['op_CC1_1_tc', 'f', 0.01, 'crystalThickness'],
    ['op_CC1_1_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_CC1_1_nvx', 'f', -0.94047583956, 'nvx'],
    ['op_CC1_1_nvy', 'f', 6.39043618484e-09, 'nvy'],
    ['op_CC1_1_nvz', 'f', -0.339860552586, 'nvz'],
    ['op_CC1_1_tvx', 'f', -0.339860552586, 'tvx'],
    ['op_CC1_1_tvy', 'f', 2.30931734946e-09, 'tvy'],
    ['op_CC1_1_ang', 'f', 0.346768620299, 'grazingAngle'],
    ['op_CC1_1_amp_coef', 'f', 1.0, 'heightAmplification'],

    # CC1_1_CC1_2: drift
    ['op_CC1_1_CC1_2_L', 'f', 0.0127, 'length'],

    # CC1_2: crystal
    ['op_CC1_2_hfn', 's', '', 'heightProfileFile'],
    ['op_CC1_2_dim', 's', 'x', 'orientation'],
    ['op_CC1_2_d_sp', 'f', 1.9201374688, 'dSpacing'],
    ['op_CC1_2_psi0r', 'f', -1.08291732223e-05, 'psi0r'],
    ['op_CC1_2_psi0i', 'f', 1.82929605452e-07, 'psi0i'],
    ['op_CC1_2_psiHr', 'f', -6.58335380584e-06, 'psiHr'],
    ['op_CC1_2_psiHi', 'f', 1.76831595902e-07, 'psiHi'],
    ['op_CC1_2_psiHBr', 'f', -6.58335380584e-06, 'psiHBr'],
    ['op_CC1_2_psiHBi', 'f', 1.76831595902e-07, 'psiHBi'],
    ['op_CC1_2_tc', 'f', 0.01, 'crystalThickness'],
    ['op_CC1_2_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_CC1_2_nvx', 'f', 0.94047583956, 'nvx'],
    ['op_CC1_2_nvy', 'f', 6.39043618484e-09, 'nvy'],
    ['op_CC1_2_nvz', 'f', -0.339860552586, 'nvz'],
    ['op_CC1_2_tvx', 'f', 0.339860552586, 'tvx'],
    ['op_CC1_2_tvy', 'f', 2.30931734946e-09, 'tvy'],
    ['op_CC1_2_ang', 'f', 0.346768620299, 'grazingAngle'],
    ['op_CC1_2_amp_coef', 'f', 1.0, 'heightAmplification'],

    # CC1_2_CC2_1: drift
    ['op_CC1_2_CC2_1_L', 'f', 0.6604, 'length'],
    #['op_Watchpoint_CC2_1_L', 'f', 0.6604, 'length'],

    # Use this to introduce an angle to CC branch:
    ['op_CC_ang_y', 'f', -15.4e-06, 'angle'],

    # CC2_1: crystal
    ['op_CC2_1_hfn', 's', '', 'heightProfileFile'],
    ['op_CC2_1_dim', 's', 'x', 'orientation'],
    ['op_CC2_1_d_sp', 'f', 1.9201374688, 'dSpacing'],
    ['op_CC2_1_psi0r', 'f', -1.08291732223e-05, 'psi0r'],
    ['op_CC2_1_psi0i', 'f', 1.82929605452e-07, 'psi0i'],
    ['op_CC2_1_psiHr', 'f', -6.58335380584e-06, 'psiHr'],
    ['op_CC2_1_psiHi', 'f', 1.76831595902e-07, 'psiHi'],
    ['op_CC2_1_psiHBr', 'f', -6.58335380584e-06, 'psiHBr'],
    ['op_CC2_1_psiHBi', 'f', 1.76831595902e-07, 'psiHBi'],
    ['op_CC2_1_tc', 'f', 0.01, 'crystalThickness'],
    ['op_CC2_1_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_CC2_1_nvx', 'f', 0.94047583956, 'nvx'],
    ['op_CC2_1_nvy', 'f', 6.39043618484e-09, 'nvy'],
    ['op_CC2_1_nvz', 'f', -0.339860552586, 'nvz'],
    ['op_CC2_1_tvx', 'f', 0.339860552586, 'tvx'],
    ['op_CC2_1_tvy', 'f', 2.30931734946e-09, 'tvy'],
    ['op_CC2_1_ang', 'f', 0.346768620299, 'grazingAngle'],
    ['op_CC2_1_amp_coef', 'f', 1.0, 'heightAmplification'],

    # CC2_1_CC2_2: drift
    ['op_CC2_1_CC2_2_L', 'f', 0.0127, 'length'],

    # CC2_2: crystal
    ['op_CC2_2_hfn', 's', '', 'heightProfileFile'],
    ['op_CC2_2_dim', 's', 'x', 'orientation'],
    ['op_CC2_2_d_sp', 'f', 1.9201374688, 'dSpacing'],
    ['op_CC2_2_psi0r', 'f', -1.08291732223e-05, 'psi0r'],
    ['op_CC2_2_psi0i', 'f', 1.82929605452e-07, 'psi0i'],
    ['op_CC2_2_psiHr', 'f', -6.58335380584e-06, 'psiHr'],
    ['op_CC2_2_psiHi', 'f', 1.76831595902e-07, 'psiHi'],
    ['op_CC2_2_psiHBr', 'f', -6.58335380584e-06, 'psiHBr'],
    ['op_CC2_2_psiHBi', 'f', 1.76831595902e-07, 'psiHBi'],
    ['op_CC2_2_tc', 'f', 0.01, 'crystalThickness'],
    ['op_CC2_2_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_CC2_2_nvx', 'f', -0.94047583956, 'nvx'],
    ['op_CC2_2_nvy', 'f', 6.39043618484e-09, 'nvy'],
    ['op_CC2_2_nvz', 'f', -0.339860552586, 'nvz'],
    ['op_CC2_2_tvx', 'f', -0.339860552586, 'tvx'],
    ['op_CC2_2_tvy', 'f', 2.30931734946e-09, 'tvy'],
    ['op_CC2_2_ang', 'f', 0.346768620299, 'grazingAngle'],
    ['op_CC2_2_amp_coef', 'f', 1.0, 'heightAmplification'],

    # Recomb_CC: aperture
    ['op_Recomb_CC_shape', 's', 'r', 'shape'],
    ['op_Recomb_CC_Dx', 'f', 0.01, 'horizontalSize'],
    ['op_Recomb_CC_Dy', 'f', 0.004, 'verticalSize'],
    ['op_Recomb_CC_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Recomb_CC_y', 'f', 0.002, 'verticalOffset'],

    ##VCC Branch
    # Split_VCC: aperture
    ['op_Split_VCC_shape', 's', 'r', 'shape'],
    ['op_Split_VCC_Dx', 'f', 0.01, 'horizontalSize'],
    ['op_Split_VCC_Dy', 'f', 0.004, 'verticalSize'],
    ['op_Split_VCC_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Split_VCC_y', 'f', -0.002, 'verticalOffset'],

    # Split_VCC_VCC1_1: drift
    ['op_Split_VCC_VCC1_1_L', 'f', 0.2, 'length'],

    # VCC1_1: crystal
    ['op_VCC1_1_hfn', 's', '', 'heightProfileFile'],
    ['op_VCC1_1_dim', 's', 'x', 'orientation'],
    ['op_VCC1_1_d_sp', 'f', 1.9201374688, 'dSpacing'],
    ['op_VCC1_1_psi0r', 'f', -1.08291732223e-05, 'psi0r'],
    ['op_VCC1_1_psi0i', 'f', 1.82929605452e-07, 'psi0i'],
    ['op_VCC1_1_psiHr', 'f', -6.58335380584e-06, 'psiHr'],
    ['op_VCC1_1_psiHi', 'f', 1.76831595902e-07, 'psiHi'],
    ['op_VCC1_1_psiHBr', 'f', -6.58335380584e-06, 'psiHBr'],
    ['op_VCC1_1_psiHBi', 'f', 1.76831595902e-07, 'psiHBi'],
    ['op_VCC1_1_tc', 'f', 0.01, 'crystalThickness'],
    ['op_VCC1_1_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_VCC1_1_nvx', 'f', 0.94047583956, 'nvx'],
    ['op_VCC1_1_nvy', 'f', 6.39043618484e-09, 'nvy'],
    ['op_VCC1_1_nvz', 'f', -0.339860552586, 'nvz'],
    ['op_VCC1_1_tvx', 'f', 0.339860552586, 'tvx'],
    ['op_VCC1_1_tvy', 'f', 2.30931734946e-09, 'tvy'],
    ['op_VCC1_1_ang', 'f', 0.346768620299, 'grazingAngle'],
    ['op_VCC1_1_amp_coef', 'f', 1.0, 'heightAmplification'],

    # VCC1_1_VCC1_2: drift
    ['op_VCC1_1_VCC1_2_L', 'f', 0.0127, 'length'],

    # VCC1_2: crystal
    ['op_VCC1_2_hfn', 's', '', 'heightProfileFile'],
    ['op_VCC1_2_dim', 's', 'x', 'orientation'],
    ['op_VCC1_2_d_sp', 'f', 1.9201374688, 'dSpacing'],
    ['op_VCC1_2_psi0r', 'f', -1.08291732223e-05, 'psi0r'],
    ['op_VCC1_2_psi0i', 'f', 1.82929605452e-07, 'psi0i'],
    ['op_VCC1_2_psiHr', 'f', -6.58335380584e-06, 'psiHr'],
    ['op_VCC1_2_psiHi', 'f', 1.76831595902e-07, 'psiHi'],
    ['op_VCC1_2_psiHBr', 'f', -6.58335380584e-06, 'psiHBr'],
    ['op_VCC1_2_psiHBi', 'f', 1.76831595902e-07, 'psiHBi'],
    ['op_VCC1_2_tc', 'f', 0.01, 'crystalThickness'],
    ['op_VCC1_2_ang_as', 'f', 0.0872664625997, 'asymmetryAngle'], #-0.0872664625997, 'asymmetryAngle'],
    ['op_VCC1_2_nvx', 'f', -0.907277635018, 'nvx'],
    ['op_VCC1_2_nvy', 'f', 6.16485781413e-09, 'nvy'],
    ['op_VCC1_2_nvz', 'f', -0.420532154534, 'nvz'],
    ['op_VCC1_2_tvx', 'f', -0.420532154534, 'tvx'],
    ['op_VCC1_2_tvy', 'f', 2.85747255185e-09, 'tvy'],
    ['op_VCC1_2_ang', 'f', 0.434031780295, 'grazingAngle'],
    ['op_VCC1_2_amp_coef', 'f', 1.0, 'heightAmplification'],

    # VCC1_2_VCC2_1: drift
    ['op_VCC1_2_VCC2_1_L', 'f', 0.2604, 'length'],
    
    # Use this to set delay between pulses:
    ['op_VCC1_2_VCC2_1_dt', 'f', 0.e-15, 'aux. time delay in VCC branch [s]'],
    #['op_VCC1_2_VCC2_1_dt', 'f', 10.e-15, 'aux. time delay in VCC branch [s]'],
    #['op_VCC1_2_VCC2_1_dt', 'f', 100.e-15, 'aux. time delay in VCC branch [s]'],
    
    ['op_VCC_shift_x', 'f', 0.e-03, 'aux. horizontal shift of beam in VCC branch [m]'],
    ['op_VCC_shift_y', 'f', 0.e-03, 'aux. vertical shift of beam in VCC branch [m]'],

    # VCC2_1: crystal
    ['op_VCC2_1_hfn', 's', '', 'heightProfileFile'],
    ['op_VCC2_1_dim', 's', 'x', 'orientation'],
    ['op_VCC2_1_d_sp', 'f', 1.9201374688, 'dSpacing'],
    ['op_VCC2_1_psi0r', 'f', -1.08291732223e-05, 'psi0r'],
    ['op_VCC2_1_psi0i', 'f', 1.82929605452e-07, 'psi0i'],
    ['op_VCC2_1_psiHr', 'f', -6.58335380584e-06, 'psiHr'],
    ['op_VCC2_1_psiHi', 'f', 1.76831595902e-07, 'psiHi'],
    ['op_VCC2_1_psiHBr', 'f', -6.58335380584e-06, 'psiHBr'],
    ['op_VCC2_1_psiHBi', 'f', 1.76831595902e-07, 'psiHBi'],
    ['op_VCC2_1_tc', 'f', 0.01, 'crystalThickness'],
    ['op_VCC2_1_ang_as', 'f', 0.0, 'asymmetryAngle'],
    ['op_VCC2_1_nvx', 'f', -0.94047583956, 'nvx'],
    ['op_VCC2_1_nvy', 'f', 6.39043618484e-09, 'nvy'],
    ['op_VCC2_1_nvz', 'f', -0.339860552586, 'nvz'],
    ['op_VCC2_1_tvx', 'f', -0.339860552586, 'tvx'],
    ['op_VCC2_1_tvy', 'f', 2.30931734946e-09, 'tvy'],
    ['op_VCC2_1_ang', 'f', 0.346768620299, 'grazingAngle'],
    ['op_VCC2_1_amp_coef', 'f', 1.0, 'heightAmplification'],

    # VCC2_1_VCC2_2: drift
    ['op_VCC2_1_VCC2_2_L', 'f', 0.0127, 'length'],

    # VCC2_2: crystal
    ['op_VCC2_2_hfn', 's', '', 'heightProfileFile'],
    ['op_VCC2_2_dim', 's', 'x', 'orientation'],
    ['op_VCC2_2_d_sp', 'f', 1.9201374688, 'dSpacing'],
    ['op_VCC2_2_psi0r', 'f', -1.08291732223e-05, 'psi0r'],
    ['op_VCC2_2_psi0i', 'f', 1.82929605452e-07, 'psi0i'],
    ['op_VCC2_2_psiHr', 'f', -6.58335380584e-06, 'psiHr'],
    ['op_VCC2_2_psiHi', 'f', 1.76831595902e-07, 'psiHi'],
    ['op_VCC2_2_psiHBr', 'f', -6.58335380584e-06, 'psiHBr'],
    ['op_VCC2_2_psiHBi', 'f', 1.76831595902e-07, 'psiHBi'],
    ['op_VCC2_2_tc', 'f', 0.01, 'crystalThickness'],
    ['op_VCC2_2_ang_as', 'f', 0.0872664625997, 'asymmetryAngle'],
    ['op_VCC2_2_nvx', 'f', 0.966516455031, 'nvx'],
    ['op_VCC2_2_nvy', 'f', 6.56737947715e-09, 'nvy'],
    ['op_VCC2_2_nvz', 'f', -0.256604641723, 'nvz'],
    ['op_VCC2_2_tvx', 'f', 0.256604641723, 'tvx'],
    ['op_VCC2_2_tvy', 'f', 1.74360203494e-09, 'tvy'],
    ['op_VCC2_2_ang', 'f', 0.259507570412, 'grazingAngle'],
    ['op_VCC2_2_amp_coef', 'f', 1.0, 'heightAmplification'],

    # VCC2_2_Recomb_VCC: drift
    ['op_VCC2_2_Recomb_VCC_L', 'f', 0.2, 'length'],

    # Recomb_VCC: aperture
    ['op_Recomb_VCC_shape', 's', 'r', 'shape'],
    ['op_Recomb_VCC_Dx', 'f', 0.01, 'horizontalSize'],
    ['op_Recomb_VCC_Dy', 'f', 0.004, 'verticalSize'],
    ['op_Recomb_VCC_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Recomb_VCC_y', 'f', -0.002, 'verticalOffset'],

    ##Common Part After Recomb.
    # Recomb_CRL: drift
    ['op_Recomb_CRL_L', 'f', 5.6, 'length'],

    # CRL1: crl
    ['op_CRL1_foc_plane', 'f', 3, 'focalPlane'],
    ['op_CRL1_delta', 'f', 3.775973e-06, 'refractiveIndex'],
    ['op_CRL1_atten_len', 'f', 0.00843, 'attenuationLength'],
    ['op_CRL1_shape', 'f', 1, 'shape'],
    ['op_CRL1_apert_h', 'f', 0.002, 'horizontalApertureSize'],
    ['op_CRL1_apert_v', 'f', 0.002, 'verticalApertureSize'],
    ['op_CRL1_r_min', 'f', 5e-05, 'tipRadius'],
    ['op_CRL1_wall_thick', 'f', 8e-05, 'tipWallThickness'],
    ['op_CRL1_x', 'f', 0.0, 'horizontalOffset'],
    ['op_CRL1_y', 'f', 0.0, 'verticalOffset'],
    ['op_CRL1_n', 'i', 7, 'numberOfLenses'],
    # CRL2: crl
    ['op_CRL2_foc_plane', 'f', 3, 'focalPlane'],
    ['op_CRL2_delta', 'f', 3.775973e-06, 'refractiveIndex'],
    ['op_CRL2_atten_len', 'f', 0.00843, 'attenuationLength'],
    ['op_CRL2_shape', 'f', 1, 'shape'],
    ['op_CRL2_apert_h', 'f', 0.002, 'horizontalApertureSize'],
    ['op_CRL2_apert_v', 'f', 0.002, 'verticalApertureSize'],
    ['op_CRL2_r_min', 'f', 0.0001, 'tipRadius'],
    ['op_CRL2_wall_thick', 'f', 8e-05, 'tipWallThickness'],
    ['op_CRL2_x', 'f', 0.0, 'horizontalOffset'],
    ['op_CRL2_y', 'f', 0.0, 'verticalOffset'],
    ['op_CRL2_n', 'i', 1, 'numberOfLenses'],
    # CRL3: crl
    ['op_CRL3_foc_plane', 'f', 3, 'focalPlane'],
    ['op_CRL3_delta', 'f', 3.775973e-06, 'refractiveIndex'],
    ['op_CRL3_atten_len', 'f', 0.00843, 'attenuationLength'],
    ['op_CRL3_shape', 'f', 1, 'shape'],
    ['op_CRL3_apert_h', 'f', 0.002, 'horizontalApertureSize'],
    ['op_CRL3_apert_v', 'f', 0.002, 'verticalApertureSize'],
    ['op_CRL3_r_min', 'f', 0.0002, 'tipRadius'],
    ['op_CRL3_wall_thick', 'f', 8e-05, 'tipWallThickness'],
    ['op_CRL3_x', 'f', 0.0, 'horizontalOffset'],
    ['op_CRL3_y', 'f', 0.0, 'verticalOffset'],
    ['op_CRL3_n', 'i', 1, 'numberOfLenses'],

    # CRL_Monitor: drift
    ['op_CRL_Monitor_L', 'f', 0.85607 - 0.002, 'length'],

#---Propagation parameters
    ##Common Part Before Split
    ['op_Aperture_pp', 'f',              [0, 0, 1.0, 0, 0, 1.0, 3.0, 1.0, 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Aperture'],
    #['op_Aperture_pp', 'f',              [0, 0, 1.0, 0, 0, 1.0, 3.0, 1.0, 24.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Aperture'],
    #['op_Aperture_pp', 'f',              [0, 0, 1.0, 0, 0, 1.0, 3.0, 1.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Aperture'],
    #['op_Aperture_pp', 'f',              [0, 0, 1.0, 0, 0, 1.0, 3.0, 1.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Aperture'],
    #['op_Aperture_pp', 'f',              [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Aperture'],
    ['op_Aperture_Split_pp', 'f',        [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Aperture_Split_VCC'],

    ##CC Branch
    ['op_Split_CC_pp', 'f',              [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Split_CC'],
    #['op_Split_CC_pp', 'f',              [0, 0, 1.0, 0, 0, 1.0, 3.0, 1.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Split_CC'],
    ['op_CC1_1_pp', 'f',                 [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CC1_1'],
    ['op_CC1_1_CC1_2_pp', 'f',           [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CC1_1_CC1_2'],
    ['op_CC1_2_pp', 'f',                 [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CC1_2'],
    ['op_CC1_2_CC2_1_pp', 'f',           [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CC1_2_CC2_1'],

    ['op_CC_ang_pp', 'f',                [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CC_ang'],

    ['op_CC2_1_pp', 'f',                 [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CC2_1'],
    ['op_CC2_1_CC2_2_pp', 'f',           [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CC2_1_CC2_2'],
    ['op_CC2_2_pp', 'f',                 [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CC2_2'],
    ['op_Recomb_CC_pp', 'f',             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Recomb_CC'],
    
    ##VCC Branch
    ['op_Split_VCC_pp', 'f',             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Split_VCC'],
    #['op_Split_VCC_pp', 'f',             [0, 0, 1.0, 0, 0, 1.0, 3.0, 1.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Split_VCC'],
    ['op_Split_VCC_VCC1_1_pp', 'f',      [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Split_VCC_VCC1_1'],
    ['op_VCC1_1_pp', 'f',                [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VCC1_1'],
    ['op_VCC1_1_VCC1_2_pp', 'f',         [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VCC1_1_VCC1_2'],
    ['op_VCC1_2_pp', 'f',                [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VCC1_2'],
    ['op_VCC1_2_VCC2_1_pp', 'f',         [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VCC1_2_VCC2_1'],
    ['op_VCC1_2_VCC2_1_dt_pp', 'f',      [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VCC1_2_VCC2_1_dt'],
    
    ['op_VCC_shift_pp', 'f',             [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'op_VCC_shift'],
    
    ['op_VCC2_1_pp', 'f',                [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VCC2_1'],
    ['op_VCC2_1_VCC2_2_pp', 'f',         [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VCC2_1_VCC2_2'],
    ['op_VCC2_2_pp', 'f',                [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VCC2_2'],
    ['op_VCC2_2_Recomb_VCC_pp', 'f',     [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VCC2_2_Recomb_VCC'],
    ['op_Recomb_VCC_pp', 'f',            [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Recomb_VCC'],
    
    ##Common Part After Recomb.
    ['op_Recomb_CRL_pp', 'f',            [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Recomb_CRL'],
    ['op_CRL1_pp', 'f',                  [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL1'],
    ['op_CRL2_pp', 'f',                  [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL2'],
    ['op_CRL3_pp', 'f',                  [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL3'],
    #['op_CRL_pp', 'f',                   [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL'],
    ['op_CRL_Monitor_pp', 'f',           [0, 0, 1.0, 4, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL_Monitor'],
    #['op_CRL_Monitor_pp', 'f',           [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL_Monitor'],
    ['op_fin_pp', 'f',                   [0, 0, 1.0, 0, 0, 0.3, 1.0, 0.3, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'final post-propagation (resize) parameters'],

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

#------------------------------------------------------------------------------
def save_pulse_data(_wfr, _ec, _xc, _yc, _fnPrefix, _do_integ=True, _do_cuts=True, _do_plot=False, _do_multi_en=False, _nmDir=''):

    sufCA = ''
    labelPosAng = ''
    unitPosAng = ''
    unitPer = ''
    mesCA = ''
    if(_wfr.presCA == 0):
        sufCA = 'c'
        labelPosAng = 'Position'
        unitPosAng = 'm'
        unitPer = 'mm'
        mesCA = 'coordinate'
    else:
        sufCA = 'a'
        labelPosAng = 'Angle'
        unitPosAng = 'rad'
        unitPer = 'mrad'
        mesCA = 'angular'

    sufTE = ''
    labelArg1 = ''
    unitArg1 = ''
    #prefSpec = ''
    onAxSpec = ''
    integSpec = ''
    onAxFN = ''
    integFN = ''
    unitPow = ''
    mesFT = ''
    if(_wfr.presFT == 0):
        sufTE = 'e' #'_e_'
        labelArg1 = 'Photon Energy'
        unitArg1 = 'eV'
        unitPow = 'J/eV'
        #prefSpec = 'Spectral '
        onAxSpec = 'Spectral Fluence'
        integSpec = 'Spectral Energy'
        onAxFN = 'spec_fluence'
        integFN = 'spec_en'       
        mesFT = 'frequency'
    else:
        sufTE = 't' #'_t_'
        labelArg1 = 'Time'
        unitArg1 = 's'
        unitPow = 'W'
        onAxSpec = 'Power Density'
        integSpec = 'Power'
        onAxFN = 'pow_den'
        integFN = 'pow'       
        mesFT = 'time'

    print('Saving ' + mesFT + '-' + mesCA + '-domain radiation characteristics to files (for viewing/debugging): ', end='')

    t0 = time()

    fnPrefixTot = _fnPrefix + '_' + sufTE + sufCA + '_'

    intType = 0 #5 #0
    #0- "Single-Electron" Intensity; 
    #5- Re(E): Real part of Single-Electron Electric Field;
    #6- Im(E): Imaginary part of Single-Electron Electric Field;
    #7- "Single-Electron" Intensity, integrated over Time or Photon Energy (i.e. Fluence);

    meshP = deepcopy(_wfr.mesh)
    meshPDvsXY = deepcopy(meshP); meshPDvsXY.ne = 1
    meshPDvsXY.eStart = _ec; meshPDvsXY.eFin = _ec

    arPDvsT = None; arPvsT = None
    if(meshP.ne > 1):

        arPDvsT = array('f', [0]*meshP.ne)
        meshPvsT = deepcopy(meshP); meshPvsT.nx = 1; meshPvsT.ny = 1
        meshPvsT.xStart = _xc; meshPvsT.xFin = _xc
        meshPvsT.yStart = _yc; meshPvsT.yFin = _yc

        if(_do_cuts):
            #Saving On-Axis (Spectral) Power Density vs Time (or Photon Energy)
            srwlpy.CalcIntFromElecField(arPDvsT, _wfr, 6, intType, 0, _ec, _xc, _yc)
            srwlib.srwl_uti_save_intens_ascii(arPDvsT, meshPvsT, os.path.join(os.getcwd(), _nmDir, fnPrefixTot + onAxFN + '.dat'),
                                              _arLabels=[labelArg1, 'Horizontal ' + labelPosAng, 'Vertical ' + labelPosAng, onAxSpec],
                                              _arUnits=[unitArg1, unitPosAng, unitPosAng, unitPow + '/' + unitPer + '^2'])
            if(_do_plot):
                uti_plot1d(arPDvsT, [meshPvsT.eStart, meshPvsT.eFin, meshPvsT.ne],
                           labels=[labelArg1, onAxSpec, onAxSpec + ' (' + _fnPrefix + ')'],
                           units=[unitArg1, unitPow + '/' + unitPer + '^2'])

        if(_do_integ):
            #Saving (Spectral) Power vs Time (or Photon Energy)
            arPvsT = array('f', [0]*meshP.ne)
            srwlpy.CalcIntFromElecField(arPvsT, _wfr, 6, 2, 0, _ec, _xc, _yc)
            srwlib.srwl_uti_save_intens_ascii(arPvsT, meshPvsT, os.path.join(os.getcwd(), _nmDir, fnPrefixTot + integFN + '.dat'),
                                              _arLabels=[labelArg1, 'Horizontal ' + labelPosAng, 'Vertical ' + labelPosAng, integSpec],
                                              _arUnits=[unitArg1, unitPosAng, unitPosAng, unitPow])
##        #Saving Input Power Density vs T (or E) & X & Y
##        arPvsEXY = array('f', [0]*meshP.ne*meshP.nx*meshP.ny)
##        srwlpy.CalcIntFromElecField(arPvsEXY, _wfr, 6, intType, 6, _ec, _xc, _yc)
##        srwlib.srwl_uti_save_intens_ascii(arPvsEXY, meshP, os.path.join(os.getcwd(), _nmDir, fnPrefixTot + onAxFN + '_3d.dat'),
##                                          _arLabels=[labelArg1, 'Horizontal ' + labelPosAng, 'Vertical ' + labelPosAng, onAxSpec],
##                                          _arUnits=[unitArg1, unitPosAng, unitPosAng, unitPow + '/' + unitPer + '^2'])
            #Saving Fluence (i.e. Power Density Integrated over Time or Spectral Fluence Integrated over Photon Energy) vs X&Y
            arFvsXY = array('f', [0]*meshP.nx*meshP.ny)
            srwlpy.CalcIntFromElecField(arFvsXY, _wfr, 6, 7, 3, _ec, _xc, _yc)
            srwlib.srwl_uti_save_intens_ascii(arFvsXY, meshPDvsXY, os.path.join(os.getcwd(), _nmDir, fnPrefixTot + 'fluence_xy.dat'),
                                              _arLabels=[labelArg1, 'Horizontal ' + labelPosAng, 'Vertical ' + labelPosAng, 'Fluence'],
                                              _arUnits=[unitArg1, unitPosAng, unitPosAng, 'J/' + unitPer + '^2'])
            if(_do_plot):
                uti_plot1d(arPvsT, [meshPvsT.eStart, meshPvsT.eFin, meshPvsT.ne],
                           labels=[labelArg1, integSpec, integSpec + ' (' + _fnPrefix + ')'],
                           units=[unitArg1, unitPow])
                uti_plot2d1d(arFvsXY,
                             [meshPDvsXY.xStart, meshPDvsXY.xFin, meshPDvsXY.nx],
                             [meshPDvsXY.yStart, meshPDvsXY.yFin, meshPDvsXY.ny], x=0, y=0,
                             labels=['Horizontal ' + labelPosAng, 'Vertical ' + labelPosAng, 'Fluence (' + _fnPrefix + ')'],
                             units=[unitPosAng, unitPosAng, 'J/' + unitPer + '^2'])

        #Saving (Spectral) Power Density In Pulse (Spectrum) vs X and Y at different Energies or Time moments
        eSaveStepFact = 1 #5 #10 #to steer
        eStep = (meshP.eFin - meshP.eStart)/(meshP.ne - 1)
        eStepSave = eStep*eSaveStepFact
        arPDvsXY = array('f', [0]*meshP.nx*meshP.ny)
        ePh = meshP.eStart
        if(_do_multi_en and _do_cuts):
        #if(_do_multi_en == True):
            for ie in range(meshP.ne):
                srwlpy.CalcIntFromElecField(arPDvsXY, _wfr, 6, intType, 3, ePh, 0, 0)
                meshPDvsXY.eStart = ePh
                meshPDvsXY.eFin = ePh
                #srwl_uti_save_intens_ascii(arPDvsXY, meshPDvsXY, os.path.join(os.getcwd(), _nmDir, fnPrefixTot + 'pow_den_xy_' + repr(round(ePh, 5)) + unitArg1 + '.dat'),
                srwlib.srwl_uti_save_intens_ascii(arPDvsXY, meshPDvsXY, os.path.join(os.getcwd(), _nmDir, fnPrefixTot + 'pow_den_xy_' + repr(ePh) + unitArg1 + '.dat'),
                                                  _arLabels=[labelArg1, 'Horizontal ' + labelPosAng, 'Vertical ' + labelPosAng, onAxSpec],
                                                  _arUnits=[unitArg1, unitPosAng, unitPosAng, unitPow + '/' + unitPer + '^2'])
                ePh += eStepSave
                if(ePh > meshP.eFin): break

    if(_do_cuts):
        #Saving Power Density (/ Spectral Fluence) In Pulse (/ Spectrum) Center vs X and Y
        arPDvsXY = array('f', [0]*meshP.nx*meshP.ny)
        srwlpy.CalcIntFromElecField(arPDvsXY, _wfr, 6, intType, 3, _ec, 0, 0)
        srwlib.srwl_uti_save_intens_ascii(arPDvsXY, meshPDvsXY, os.path.join(os.getcwd(), _nmDir, fnPrefixTot + onAxFN + '_xy.dat'),
                                          _arLabels=[labelArg1, 'Horizontal ' + labelPosAng, 'Vertical ' + labelPosAng, onAxSpec],
                                          _arUnits=[unitArg1, unitPosAng, unitPosAng, unitPow + '/' + unitPer + '^2'])
        if(_do_plot):
            uti_plot2d1d(arPDvsXY,
                         [meshPDvsXY.xStart, meshPDvsXY.xFin, meshPDvsXY.nx],
                         [meshPDvsXY.yStart, meshPDvsXY.yFin, meshPDvsXY.ny], x=0, y=0,
                         labels=['Horizontal ' + labelPosAng, 'Vertical ' + labelPosAng, onAxSpec + ' (' + _fnPrefix + ')'],
                         units=[unitPosAng, unitPosAng, unitPow + '/' + unitPer + '^2'])

    print('done in', round(time() - t0), 's')

#------------------------------------------------------------------------------
def main(_delay=0.e-15, _use_crl=True, _do_integ=True, _do_cuts=True):
    #_do_integ #produce integrated pulse characteristics
    #_do_cuts #produce cuts of 3D pulse characteristics
    
    v = srwl_bl.srwl_uti_parse_options(varParam, use_sys_argv=True)

    #Delay of pulse in VCC branch:
    v.op_VCC1_2_VCC2_1_dt = _delay

    #Necessary for initial wavefront calculation without ropagation
    v.si = True
    
    #srwl_bl.SRWLBeamline(_name=v.name, _mag_approx=mag).calc_all(v, op)
    srwl_bl.SRWLBeamline(_name=v.name).calc_all(v)
    wfr = v.w_res

    #Coordinates of central point for resulting distribution cuts: 
    tc = 0.5*(v.w_e + v.w_ef)
    xc = v.w_x
    yc = v.w_y

    #Initial radiaiton in Time-Coordinate domain
    save_pulse_data(wfr, tc, xc, yc, _fnPrefix='init', _do_integ=_do_integ, _do_cuts=_do_cuts, _do_plot=True, _do_multi_en=False, _nmDir=v.fdir)

    print('Resizing in Time domain: ', end='')
    t0 = time();
    srwlpy.ResizeElecField(wfr, 't', [0, 25., 1])
    #srwlpy.ResizeElecField(wfr, 't', [0, 50., 1])
    print('done in', round(time() - t0, 3), 's')
    
    print('Switching from Time to Frequency domain: ', end='')
    t0 = time();
    srwlpy.SetRepresElecField(wfr, 'f');
    print('done in', round(time() - t0, 3), 's')

    print('Resizing in Frequency domain: ', end='')
    t0 = time();
    srwlpy.ResizeElecField(wfr, 'f', [0, 0.07, 1])
    print('done in', round(time() - t0, 3), 's')
    
    ec = 0.5*(wfr.mesh.eStart + wfr.mesh.eFin) #- 0.04
    
    save_pulse_data(wfr, ec, xc, yc, _fnPrefix='init', _do_integ=_do_integ, _do_cuts=_do_cuts, _do_plot=True, _do_multi_en=False, _nmDir=v.fdir)

    print('Propagating pulse to Split Point: ', end='')
    t0 = time()
    blBeforeSplit = set_optics_before_split(v)
    srwlpy.PropagElecField(wfr, blBeforeSplit)
    print('done in', round(time() - t0, 3), 's')
    
    opOrientData = blBeforeSplit.get_orient_lab_fr(_e=v.gbm_ave, _r=v.op_r) 
    uti_io.write_ascii_data_rows(_file_path=os.path.join(v.fdir, 'orient_before_split.dat'), _rows=opOrientData, _str_sep='\t', 
                                 _str_head='#Types of optical elements and Cartesian coordinates of their center positions and base vectors (t, s, n) in the Lab frame')

    save_pulse_data(wfr, ec, xc, yc, _fnPrefix='at_split', _do_integ=_do_integ, _do_cuts=_do_cuts, _do_plot=True, _do_multi_en=False, _nmDir=v.fdir)

    wfrCC = wfr
    wfrVCC = deepcopy(wfr)
    #Test:
    #wfrVCC = wfr
    
    print('Propagating half-pulse through CC Branch to Recombination Point: ', end='')
    t0 = time()
    blCC = set_optics_cc(v)
    srwlpy.PropagElecField(wfrCC, blCC)
    print('done in', round(time() - t0, 3), 's')
    
    opOrientData = blCC.get_orient_lab_fr(_e=v.gbm_ave, _r=(v.op_r + v.op_Aperture_Split_L) ) 
    uti_io.write_ascii_data_rows(_file_path=os.path.join(v.fdir, 'orient_cc.dat'), _rows=opOrientData, _str_sep='\t', 
                                 _str_head='#Types of optical elements and Cartesian coordinates of their center positions and base vectors (t, s, n) in the Lab frame')

    save_pulse_data(wfrCC, ec, xc, yc, _fnPrefix='at_recomb_cc', _do_integ=_do_integ, _do_cuts=_do_cuts, _do_plot=True, _do_multi_en=False, _nmDir=v.fdir)

    print('Propagating half-pulse through VCC Branch to Recombination Point: ', end='')
    t0 = time()
    blVCC = set_optics_vcc(v)
    srwlpy.PropagElecField(wfrVCC, blVCC)
    print('done in', round(time() - t0, 3), 's')
    
    opOrientData = blVCC.get_orient_lab_fr(_e=v.gbm_ave, _r=(v.op_r + v.op_Aperture_Split_L) ) 
    uti_io.write_ascii_data_rows(_file_path=os.path.join(v.fdir, 'orient_vcc.dat'), _rows=opOrientData, _str_sep='\t', 
                                 _str_head='#Types of optical elements and Cartesian coordinates of their center positions and base vectors (t, s, n) in the Lab frame')

    save_pulse_data(wfrVCC, ec, xc, yc, _fnPrefix='at_recomb_vcc', _do_integ=_do_integ, _do_cuts=_do_cuts, _do_plot=True, _do_multi_en=False, _nmDir=v.fdir)

#    print('CC Mesh:')
#    print(wfr.mesh.ne, wfr.mesh.eStart, wfr.mesh.eFin)
#    print(wfr.mesh.nx, wfr.mesh.xStart, wfr.mesh.xFin)
#    print(wfr.mesh.ny, wfr.mesh.yStart, wfr.mesh.yFin)
#    print('VCC Mesh:')
#    print(wfrVCC.mesh.ne, wfrVCC.mesh.eStart, wfrVCC.mesh.eFin)
#    print(wfrVCC.mesh.nx, wfrVCC.mesh.xStart, wfrVCC.mesh.xFin)
#    print(wfrVCC.mesh.ny, wfrVCC.mesh.yStart, wfrVCC.mesh.yFin)

    print('Recombining Half-Pulses: ', end='')
    t0 = time()
    wfr.addE(wfrVCC)
    del wfrVCC
    print('done in', round(time() - t0, 3), 's')

    save_pulse_data(wfr, ec, xc, yc, _fnPrefix='after_recomb', _do_integ=_do_integ, _do_cuts=_do_cuts, _do_plot=True, _do_multi_en=False, _nmDir=v.fdir)

    print('Propagating Recombined Pulse to Monitor: ', end='')
    t0 = time()
    srwlpy.PropagElecField(wfr, set_optics_after_recomb(v, _use_crl=_use_crl))
    print('done in', round(time() - t0, 3), 's')

    save_pulse_data(wfr, ec, xc, yc, _fnPrefix='at_monitor', _do_integ=_do_integ, _do_cuts=_do_cuts, _do_plot=True, _do_multi_en=False, _nmDir=v.fdir)

    print('Resizing in Frequency domain: ', end='')
    t0 = time();
    srwlpy.ResizeElecField(wfr, 'f', [0, 3., 1.])
    print('done in', round(time() - t0, 3), 's')

    print('Switching from Frequency to Time domain: ', end='')
    t0 = time();
    srwlpy.SetRepresElecField(wfr, 't')
    print('done in', round(time() - t0, 3), 's')

    save_pulse_data(wfr, tc, xc, yc, _fnPrefix='at_monitor', _do_integ=_do_integ, _do_cuts=_do_cuts, _do_plot=True, _do_multi_en=False, _nmDir=v.fdir)
    
    uti_plot_show()

#------------------------------------------------------------------------------
if __name__ == '__main__':
    main(_delay=0.e-15, _use_crl=False, _do_integ=True, _do_cuts=False)
