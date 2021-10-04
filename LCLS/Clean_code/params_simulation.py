from functions_calc import *

''' source parameters '''
N = 1024           # number of sampling points
E0 = 9000          # photon energy [eV]
sigma_x = 23e-6    # RMS source size in x
sigma_y = 23e-6    # RMS source size in y
z_source = 650     # source position
beam_params = {
    'photonEnergy': E0,
    'N': N,
    'sigma_x': sigma_x,
    'sigma_y': sigma_y,
    'rangeFactor': 5,
    'scaleFactor': 10,
    'z_source': z_source
}


''' telescope mirror parameters '''
m1_p = 185.0              # telescope mirror 1 source distance [m]
m_alpha = 2.65e-3         # telescope mirror grazing angle [rad]
m2_z = 300.0              # telescope mirror 2 position (from source) [m]
d_m1_m2 = m2_z - m1_p     # distance between telescope mirrors [m]
nFWHM = 1.0               # +- nxFWHM of the beam on the first telescope mirror

m1_q, m2_p = calc_mir_param(E0, sigma_x, m1_p, nFWHM=nFWHM, m_alpha=m_alpha, d_m1_m2=d_m1_m2)


''' crystal parameters '''
HHLM_type = '2DCM'        # 2DCM or Zigzag
HHLM_offset = 20e-3       # distance between off-axis OEs and the nominal axis [m]
pair_distance=200e-3      # distance between HRM crystal 2 and 3 [m]

hkl1 = [1,1,1]; hkl2 = [4,4,0]; hkl3 = [4,4,0]    # 
hkl1 = [1,1,1]; alphaAsym1 = 9.0    # Miller planes and asymmetry angles for HHLM1,2
hkl2 = [4,4,0]; alphaAsym2 = 38.3   # ````````````````````````HHLM3,4
hkl3 = [4,4,0]; alphaAsym3 = 15.0   # ````````````````````````HRM C2, C3

crystals = np.array(['HHLM1', 'HHLM2', 'HHLM3', 'HHLM4', 'C1', 'C2']); n_crys = len(crystals)

dir_out = 'output/'; make_dir(dir_out)
dir_case = dir_out + '{}eV_{}_{}-{}_{}-{}_{}mm/'.format(E0, HHLM_type,
                                          hkl1[0]*100+hkl1[1]*10+hkl1[2], alphaAsym1,
                                          hkl2[0]*100+hkl2[1]*10+hkl2[2], alphaAsym2,
                                          HHLM_offset*1e3); make_dir(dir_case)


''' HRM mirror parameters '''
f1 = 10.0
f2 = 10.0


''' shapeErrors '''
dir_profile = '../DXS_Xtals_FEA4WFS_3x/n_fwhm_6/{}/'.format(HHLM_type)
nFWHMx_profile = 3; nFWHMy_profile = 2

if HHLM_type == 'Zigzag':
    crystal_profile_powers = np.array([
        [73.15, 16.53, 0.029, 0.027, 0.026, 0.03],
        [73.15, 16.79, 0.015, 0.015, 0.014, 0.013],
        [4.51, 50.1, 1.78, 1.69, 1.61, 1.88] ])
elif HHLM_type == '2DCM':
    crystal_profile_powers = np.array([
        [73.15, 0.027, 16.53, 0.029, 0.026, 0.03],
        [73.15, 0.015, 16.79, 0.015, 0.014, 0.013],
        [4.51, 1.69, 50.1, 1.78, 1.61, 1.88] ])


''' phase corrector '''
FOV1 = 5e-3; N1 = 512
FOV2 = 5e-3; N2 = 1024
plate_position = 'HHLM'


''' case specific parameters (variables) '''
nFWHMx = 3
nFWHMy = 1
name_size = '{}FWHMx_{}FWHMy'.format(nFWHMx, nFWHMy)

m2_p_perfect = 167.84807686607576
m2_p_bending = 2725.7090355414666

crystal_powers_case = np.array([29.8626,0.7631,13.3021,0.3261,0.0889,0.05])    # crystal power for 9keV 2DCM-111_9.0-440_38.3-440_15.0 SASE

