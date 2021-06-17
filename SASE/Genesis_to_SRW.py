from h5py import File
from genesis.writers import write_openpmd_wavefront_h5
from Module_diagnostic_functions import *

# DFL = 'stats_9.00_0.out.dfl'

# # Genesis parameters. These need to be known to populate the openPMD-wavefront metadata
# HXREn=9 #keV
# PARAM={'ncar': 251, 'dgrid': 0.0007, 'xlamds': 1240e-9/(HXREn*1000), 'zsep': 20, 'ntail':0, 'itdp':1}

# from genesis.parsers import parse_genesis_dfl
# D2 = parse_genesis_dfl(DFL, nx=PARAM['ncar'])

# print(PARAM['xlamds']*PARAM['zsep'])

# with File('wavefront.h5', 'w') as h5:
#     write_openpmd_wavefront_h5(h5, dfl=D2, param=PARAM)
# h5 = File('wavefront.h5', 'r')
# dict(h5['data']['000000']['meshes']['electricField'].attrs)
# print('SRW wavefront written')

#### To SRW
from pmd_srw import srw_wfr_from_openpmd_wavefront

with File('wavefront.h5', 'r') as h5:
    arEx, arEy, kwargs,  wfr_attrs = srw_wfr_from_openpmd_wavefront(h5['data']['000000']['meshes'],  iz_step=None)

wfr = srwlib.SRWLWfr(_arEx=arEx, _arEy=arEy, **kwargs)
wfr.__dict__.update(wfr_attrs)

##### view beam
plt.figure(figsize=(17,5))
plt.subplot(1,3,1); plot_spatial_from_wf(wfr, if_slice=0, if_log=0)
plt.subplot(1,3,2); plot_tilt_from_wf(wfr, ori='Vertical', type='sum', if_log=0)
plt.subplot(1,3,3); plot_spectrum_from_wf(wfr, if_short=1)

plt.savefig('profile.png')




# import numpy as np
# from array import array
# from copy import deepcopy
# from Module_diagnostic_functions import *

# def parse_genesis_dfl(fname, nx):
#     """
#     fname: filename
#     nx: grid size in x and y. Same as Genesis 'ncar'
    
#     Returns 3d numpy.array with indices as:
#         [x, y, z]
    
#     """
#     dat = np.fromfile(fname, dtype=np.complex).astype(np.complex)
#     npoints = dat.shape[0] 
    
#     # Determine number of slices
#     ny = nx
#     nz =  npoints / ny /nx
#     assert (nz % 1 == 0), f'Confused shape {nx} {ny} {nz}' 
#     nz = int(nz)   
#     dat = dat.reshape(nz, ny, nx)    
#     dat = np.moveaxis(dat, [0,1,2], [2,1,0]) # z, y, x to x, y, z
    
#     return dat
	
# def dfl_to_E(dfl,ncar,dgrid):
#     """
#     This function takes a dfl (from parsers.parse_genesis_dfl) ncar, and dgrid, and returns the electric field (in SI units, V/m)
#     from genesis.analysis (David's LUME)
#     """
#     #compare to intensity
#     area=(dgrid*2/(ncar-1))**2
#     return np.sqrt(2*376.7)*dfl/np.sqrt(area)
	

# #Constants from Sim
# HXREn=9 #keV
# inputfile='stats_9.00_0.out.dfl'

# ncar=251
# dgrid=0.0007
# xlamds=1240e-9/(HXREn*1000)
# zsep=20

# #Load wavefront
# Efield=parse_genesis_dfl(inputfile,ncar) #Load DFL, 3D array ordered [z,x,y]
# #Process
# pulse=dfl_to_E(Efield,ncar,dgrid) #Convert to V/m
# zs=np.asarray(range(Efield.shape[0]))*xlamds*zsep
# dxy=2*dgrid/float(ncar-1)
# xys=(np.asarray(range(ncar))-(ncar-1)/2)*dxy #ncar is odd
# xs=xys;
# ys=xys;


# #### changing format
# delta = xs[1] - xs[0]
# dz = zs[1] - zs[0]

# polarization = 'y'
# photon_energy_eV = HXREn*1e3

# nx, ny, nz = Efield.shape

# iz_start = 0; iz_end = nz; iz_step = 1
# zslice = slice(iz_start, iz_end, iz_step)

# # Input grid info
# mins = [-nx*delta/2, -ny*delta/2, -nz*dz/2]
# maxs = np.add([(nx-1)*delta, (ny-1)*delta, (nz-1)*dz], mins)

# def srwl_uti_array_alloc(_type, _n):
#     """
#     make empty array with SRW wavefront structure
#     """
#     nPartMax = 10000000 #to tune
#     if(_n <= nPartMax): return array(_type, [0]*_n)
#         #resAr = array(_type, [0]*_n)
#         #print('Array requested:', _n, 'Allocated:', len(resAr))
#         #return resAr

#     nEqualParts = int(_n/nPartMax)
#     nResid = int(_n - nEqualParts*nPartMax)
#     resAr = array(_type, [0]*nPartMax)
#     if(nEqualParts > 1):
#         auxAr = deepcopy(resAr)
#         for i in range(nEqualParts - 1): resAr.extend(auxAr)
#     if(nResid > 0):
#         auxAr = array(_type, [0]*nResid)
#         resAr.extend(auxAr)

#     #print('Array requested:', _n, 'Allocated:', len(resAr))
#     return resAr

# def np_complex_to_srw_array(cdat):
# 	"""
#     convert 3D complex field to SRW wavefront structure
#     """
# 	dat = np.moveaxis(cdat, [0,1,2], [1,0,2]).flatten() # [x,y,z] to # [y,x,z]
# 	x_np = np.ravel([np.real(dat), np.imag(dat)], order = 'F').astype(np.dtype('f'))
# 	return array(x_np.dtype.char, x_np.data.tobytes())

# if polarization == 'x':
# 	arEx = np_complex_to_srw_array(pulse[:,:,zslice])
# 	arEy = None
# else:
# 	arEx = None
# 	arEy = np_complex_to_srw_array(pulse[:,:,zslice])

# if arEx is None: arEx = srwl_uti_array_alloc('f', len(arEy))
# if arEy is None: arEy = srwl_uti_array_alloc('f', len(arEx))

# # Actual grid info
# c_light = 299792458. # m/s

# iz_end = iz_step*((iz_end-iz_start)//iz_step) + iz_start
# nz = (iz_end - iz_start)//iz_step

# zmin = iz_start*dz +mins[2]
# zmax = (iz_end-1)*dz +mins[2]

# # kwargs for SRWLWfr init
# kwargs = dict(_typeE  = 'f',
#               _eStart = 0, #zmin/c_light,
#               _eFin   = (zmax-zmin)/c_light,
#               _ne     = nz,
#               _xStart =mins[0],
#               _xFin   =maxs[0],
#               _nx     =nx,
#               _yStart =mins[1],
#               _yFin   =maxs[1],
#               _ny     =ny,
#               _zStart =0)

# # wfr. items
# presFT = 1
# wfr_attrs = dict(avgPhotEn = photon_energy_eV, presFT = presFT, unitElFld=2)

# wfr = srwlib.SRWLWfr(_arEx=arEx, _arEy=arEy, **kwargs)
# wfr.__dict__.update(wfr_attrs)

# ##### view beam
# plt.figure(figsize=(17,5))
# plt.subplot(1,3,1); plot_spatial_from_wf(wfr, if_slice=0, if_log=0)
# plt.subplot(1,3,2); plot_tilt_from_wf(wfr, ori='Vertical', type='sum', if_log=0)
# plt.subplot(1,3,3); plot_spectrum_from_wf(wfr, if_short=1)

# plt.savefig('profile.png')
