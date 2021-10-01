from functions_misc import *
from functions_beamline import *
from skimage.restoration import unwrap_phase

''' calculate mirror dimension '''
def calc_mir_param(E0, sigma_x, m1_p, nFWHM=1.0, m_alpha=2.65e-3, d_m1_m2=115.0):
    # m1_p: mirror 1 source distance [m]
    # E0: photon energy [eV]
    # sigma_x: source RMS size [m]
    # nFWHM: beam size on the 1st mirror (i.e. n=2 means beam size = 50% nominal)
    # m_alpha: telescope mirror grazing angle [rad]
    # d_m1_m2: distance between telescope mirrors [m]
    wavelength = 1239.84193e-9/E0
    Rayleigh_length = np.pi * sigma_x**2 / wavelength
    beam_width_z = 2.35 * sigma_x * np.sqrt(1+np.square(m1_p/Rayleigh_length))    # 2xFWHM

    m2_p = m_alpha/nFWHM/(m_alpha/nFWHM-beam_width_z)*d_m1_m2
    m1_q = d_m1_m2 - m2_p
    return m1_q, m2_p


''' calculate crystal dimension '''
def calc_crystal_dimension(E0, z_s=650, m_alpha=2.65e-3, m2_z=300.0,
                           m1_p=185.0, m1_q=-25.6, m2_p=141.6,
                           aperture_size=1.0,

                           HHLM_type='2DCM',
                           HHLM_offset=20e-3,
                           pair_distance=200e-3,
                           hkl1=[1,1,1], alphaAsym1=0.0,
                           hkl2=[2,2,0], alphaAsym2=0.0,
                           
                           f1=10.0, f2=10.0, slit_width=1e-1,
                           hkl3=[5,5,5], alphaAsym3=15.0,
                           beam_params=None,
                           crystals=['HHLM1','HHLM2','HHLM3','HHLM4'],
                           nFWHMx=1, nFWHMy=1,
                           stdout=None):

    # simulate parameters for nominal condition
    sys.stdout = open(os.devnull, 'w')
    mono_devices_calc, mono_beamline_calc = define_beamline_normal(
                           E0, z_s=z_s, m_alpha=m_alpha, m2_z=m2_z,
                           m1_p=m1_p, m1_q=m1_q, m2_p=m2_p,
                           aperture_size=aperture_size,

                           HHLM_type=HHLM_type,
                           HHLM_offset=HHLM_offset,
                           pair_distance=pair_distance,
                           hkl1=hkl1, alphaAsym1=alphaAsym1,
                           hkl2=hkl2, alphaAsym2=alphaAsym2,
                           
                           f1=10.0, f2=10.0, slit_width=slit_width,
                           hkl3=hkl3, alphaAsym3=alphaAsym3)

    b1_calc = beam.Beam(beam_params=beam_params)
    b2_calc = mono_beamline_calc.propagate_beamline(b1_calc)
    sys.stdout = stdout

    # get crystal properties
    beta0s = []; braggs = []; asyms = []; OE_z = []
    for device in mono_devices_calc:
        if device.name in crystals:
            beta0s.append(device.beta0)
            braggs.append(np.rad2deg(device.bragg))
            asyms.append(np.rad2deg(device.alphaAsym))
            OE_z.append(device.z)

    # calculate beam FWHM
    wx1, wy1 = mono_beamline_calc.im_after_HHLM1.beam_analysis(mono_beamline_calc.im_after_HHLM1.get_profile_x(), mono_beamline_calc.im_after_HHLM1.get_profile_y())[2:4]
    wx2, wy2 = mono_beamline_calc.im_after_HHLM2.beam_analysis(mono_beamline_calc.im_after_HHLM2.get_profile_x(), mono_beamline_calc.im_after_HHLM2.get_profile_y())[2:4]
    wx3, wy3 = mono_beamline_calc.im_after_HHLM3.beam_analysis(mono_beamline_calc.im_after_HHLM3.get_profile_x(), mono_beamline_calc.im_after_HHLM3.get_profile_y())[2:4]
    wx4, wy4 = mono_beamline_calc.im_after_HHLM4.beam_analysis(mono_beamline_calc.im_after_HHLM4.get_profile_x(), mono_beamline_calc.im_after_HHLM4.get_profile_y())[2:4]
    wx5, wy5 = mono_beamline_calc.im_after_C1.beam_analysis(mono_beamline_calc.im_after_C1.get_profile_x(), mono_beamline_calc.im_after_C1.get_profile_y())[2:4]
    wx6, wy6 = mono_beamline_calc.im_after_C2.beam_analysis(mono_beamline_calc.im_after_C2.get_profile_x(), mono_beamline_calc.im_after_C2.get_profile_y())[2:4]
    beam_size = np.array([wx1, wx2, wx3, wx4, wx5, wx6])

    # crystal dimension
    l_crystal = np.divide(beam_size*2*nFWHMx, np.sin(beta0s))
    w_crystal = np.array([wy1, wy2, wy3, wy4, wy5, wy6])*2*nFWHMy

    print('Bragg angle: {}'.format(np.round(braggs,3)))
    print('asymmetry: {}'.format(np.round(asyms,3)))
    print('beam FWHM: {} mm'.format(np.round(beam_size*1e3, 3)))
    print('footprint x: {} mm'.format(np.round(l_crystal*1e3, 3)))
    print('footprint y: {} mm'.format(np.round(w_crystal*1e3, 3)))
    print('position: {} m\n'.format(np.round(OE_z, 5)))

    return beam_size, l_crystal, w_crystal


''' calculates the width of the slit in HRM '''
def calc_slit_width(inbeam, lmbd, foc, min_width=3e-6):
    # inbeam: the beamsize (2x FWHM) after Crystal 2
    # lmbd: wavelength of photon
    # asym: asymmetry angle of Crystal 2
    # foc: Mir1/Lens1 focal distance
    # returns mono focus size in [m] at the Fourier plane
    slt = 2*np.log(2)*lmbd*foc/np.pi/inbeam
    slit_width = 2.11 * slt

    if slit_width < min_width:
        print('calculated slit width too narrow ({} um)'.format(round(slit_width*1e6, 2)))
        slit_width = min_width
    print('actual slit width: {} um\n'.format(round(slit_width*1e6, 2)))
    return slit_width


''' calculate closest crystal profile '''
def calc_crystal_option(profile_power, actual_power):
    ratio_power = profile_power/np.stack([actual_power for i in range(len(profile_power))])
    ratio_power[ratio_power<1] = 1/(ratio_power[ratio_power<1])
    options = ratio_power.argmin(axis=0)+1
    return options


''' calculate crystal shapeErrors '''
def calc_crystal_shapeErrors(dir_profile,
                             crystal_names, crystal_powers,
                             crystal_profile_powers, options,
                             l_crystal, w_crystal,
                             nFWHMx=3, nFWHMx_profile=3,
                             nFWHMy=1, nFWHMy_profile=1):
    shapeErrors = [None for i in range(len(crystal_names))]
    x_mirs = [None for i in range(len(crystal_names))]
    y_mirs = [None for i in range(len(crystal_names))]
    
    n_fig = len(crystal_names)
    plt.figure(figsize=(2+n_fig*5, 8))
    for i, crystal in enumerate(crystal_names):
        if crystal_powers[i] == 0: continue
            
        # load profiles from file
        dy, xx, zz = load_crystal_data(dir_profile, crystal, options[i],
                                       nFWHMx=nFWHMx, nFWHMx_profile=nFWHMx_profile,
                                       nFWHMy=nFWHMy, nFWHMy_profile=nFWHMy_profile)
        length_profile = xx.max() - xx.min()
        width_profile = zz.max() - zz.min()
        
        # proportionality constant due to power and crystal dimension
        ratio_power = crystal_powers[i]/crystal_profile_powers[options[i]-1, i]
        ratio_length = l_crystal[i]*1e3/ length_profile
        ratio_width = w_crystal[i]*1e3 / width_profile
        
        shapeErrors[i] = dy*1e6 * ratio_power / ratio_length / ratio_width    # [nm]
        x_mirs[i] = np.linspace(-l_crystal[i]*1e3/2, l_crystal[i]*1e3/2, 2048)
        y_mirs[i] = np.linspace(-w_crystal[i]*1e3/2, w_crystal[i]*1e3/2, 2048)
    
        plt.subplot(2,n_fig,i+1)
        plt.imshow(dy*1e6, cmap='jet',
                   extent = (xx.min(), xx.max(), zz.max(), zz.min()))
        plt.colorbar(); plt.axis('tight'); plt.title(crystal); plt.ylabel('loaded (mm)')
        plt.subplot(2,n_fig,n_fig+i+1)
        plt.imshow(shapeErrors[i], cmap='jet',
                   extent = (x_mirs[i].min(), x_mirs[i].max(), y_mirs[i].max(), y_mirs[i].min()))
        plt.colorbar(); plt.axis('tight'); plt.xlabel('x (mm)'); plt.ylabel('interpolated (mm)')
        plt.title('{}p/{}x/{}y'.format(round(ratio_power,2), round(ratio_length,2), round(ratio_width,2)))
    return shapeErrors, x_mirs, y_mirs


''' calculate phase difference '''
def calc_phase_diff(E0, z_s=650, beam_params=None,
                    m_alpha=2.65e-3, m2_z=300.0,
                    m1_p=185.0, m1_q=-25.6,
                    m2_p_perfect=141.6, m2_p_bending=2000.0,
                    aperture_size=1.0,

                    HHLM_type='2DCM',
                    HHLM_offset=20e-3,
                    pair_distance=200e-3,
                    hkl1=[1,1,1], alphaAsym1=0.0,
                    hkl2=[2,2,0], alphaAsym2=0.0,
                    l_crystal=[1e-1 for i in range(6)],
                    w_crystal=[5e-3 for i in range(6)],

                    f1=10.0, f2=10.0, slit_width=1e-1,
                    hkl3=[5,5,5], alphaAsym3=15.0,
                    shapeErrors=[None for i in range(6)],

                    FOV1=5e-3, N1=512,
                    FOV2=5e-3, N2=8192,
                    plate_position='HHLM',
                   
                    nFWHMx=3, nFWHMy=1,
                    stdout = None):
    
    # define beamlines w/ or w/o shapeError
    sys.stdout = open(os.devnull, 'w')
    _, mono_beamline_perfect = define_beamline_normal(
                            E0, z_s=z_s, m_alpha=m_alpha, m2_z=m2_z,
                            m1_p=m1_p, m1_q=m1_q, m2_p=m2_p_perfect,
                            aperture_size=aperture_size,
    
                            HHLM_type=HHLM_type,
                            HHLM_offset=HHLM_offset,
                            pair_distance=pair_distance,
                            hkl1=hkl1, alphaAsym1=alphaAsym1,
                            hkl2=hkl2, alphaAsym2=alphaAsym2,
                            l_crystal=l_crystal,
                            w_crystal=w_crystal,
                            
                            f1=f1, f2=f2, slit_width=slit_width,
                            hkl3=hkl3, alphaAsym3=alphaAsym3,
                            shapeErrors=[None for i in range(6)],
    
                            FOV1=FOV1, N1=N1,
                            FOV2=FOV2, N2=N2,
                            plate_position=plate_position)

    _, mono_beamline_bending = define_beamline_normal(
                            E0, z_s=z_s, m_alpha=m_alpha, m2_z=m2_z,
                            m1_p=m1_p, m1_q=m1_q, m2_p=m2_p_bending,
                            aperture_size=aperture_size,
    
                            HHLM_type=HHLM_type,
                            HHLM_offset=HHLM_offset,
                            pair_distance=pair_distance,
                            hkl1=hkl1, alphaAsym1=alphaAsym1,
                            hkl2=hkl2, alphaAsym2=alphaAsym2,
                            l_crystal=l_crystal,
                            w_crystal=w_crystal,
                            
                            f1=f1, f2=f2, slit_width=slit_width,
                            hkl3=hkl3, alphaAsym3=alphaAsym3,
                            shapeErrors=shapeErrors,
    
                            FOV1=FOV1, N1=N1,
                            FOV2=FOV2, N2=N2,
                            plate_position=plate_position)
    
    # propagate monochromatic results
    beam_params['photonEnergy'] = E0
    b1_perfect = beam.Beam(beam_params=beam_params)
    b2_perfect = mono_beamline_perfect.propagate_beamline(b1_perfect)
    
    beam_params['photonEnergy'] = E0
    b1_bending = beam.Beam(beam_params=beam_params)
    b2_bending = mono_beamline_bending.propagate_beamline(b1_bending)
    sys.stdout = stdout
    
    # extract dimensional info of phase plate
    ppm_perfect = mono_beamline_perfect.im_plate
    ppm_bending = mono_beamline_bending.im_plate
    
    plate_length = nFWHMx * ppm_perfect.get_x_width()
    plate_width = nFWHMy * ppm_perfect.get_y_width()
    index_x = np.where((ppm_perfect.x ** 2) < ((plate_length/2)**2))[0]
    index_y = np.where((ppm_perfect.y ** 2) < ((plate_width/2)**2))[0]
    x_plate = np.linspace(-plate_length/2, plate_length/2, index_x[-1]-index_x[0]+1) * 1e3
    y_plate = np.linspace(-plate_width/2, plate_width/2, index_y[-1]-index_y[0]+1) * 1e3
    
    # extract phase difference
    if plate_position == 'HRM':
        ppm_perfect = mono_beamline_perfect.im_plate2
        ppm_bending = mono_beamline_bending.im_plate2
    phase_perfect_x = ppm_perfect.x_phase[index_x[0]: index_x[-1]+1]
    phase_perfect_y = ppm_perfect.y_phase[index_y[0]: index_y[-1]+1]
    phase_bending_x = ppm_bending.x_phase[index_x[0]: index_x[-1]+1]
    phase_bending_y = ppm_bending.y_phase[index_y[0]: index_y[-1]+1]
    
    phase_perfect_2d = unwrap_phase(np.angle(ppm_perfect.complex_beam()[0]), wrap_around=(False, False))[index_y[0]:index_y[-1]+1, index_x[0]:index_x[-1]+1]
    phase_bending_2d = unwrap_phase(np.angle(ppm_bending.complex_beam()[0]), wrap_around=(False, False))[index_y[0]:index_y[-1]+1, index_x[0]:index_x[-1]+1]
    
    phase_perfect_2d -= np.median(phase_perfect_2d)
    phase_bending_2d -= np.median(phase_bending_2d)
    phase_diff_2d = phase_bending_2d - phase_perfect_2d
    
    # plot phase
    plt.figure(figsize=(17,5))
    plt.subplot(1,3,1)
    plt.imshow(phase_perfect_2d, cmap='jet',
               extent=[-plate_length*1e3/2, plate_length*1e3/2, plate_width*1e3/2, -plate_width*1e3/2])
    plt.colorbar(); plt.title('perfect'); plt.xlabel('x (mm)'); plt.ylabel('y (mm)'); plt.axis('tight')
    
    plt.subplot(1,3,2)
    plt.imshow(phase_bending_2d, cmap='jet',
               extent=[-plate_length*1e3/2, plate_length*1e3/2, plate_width*1e3/2, -plate_width*1e3/2])
    plt.colorbar(); plt.title('bending'); plt.xlabel('x (mm)'); plt.ylabel('y (mm)'); plt.axis('tight')
    
    plt.subplot(1,3,3)
    plt.imshow(phase_diff_2d, cmap='jet',
               extent=[-plate_length*1e3/2, plate_length*1e3/2, plate_width*1e3/2, -plate_width*1e3/2])
    plt.colorbar(); plt.title('phase difference'); plt.xlabel('x (mm)'); plt.ylabel('y (mm)'); plt.axis('tight')
    
    return phase_diff_2d, x_plate, y_plate


''' calculate phase plate thickness '''
def calc_plate_thickness(E0, phase_diff_2d, x_plate, y_plate,
                         nFWHMx=3, nFWHMx_fit=3,
                         nFWHMy=1, nFWHMy_fit=1):
    index_x = np.where((x_plate**2) < (x_plate.max()/nFWHMx*nFWHMx_fit)**2)[0]
    index_y = np.where((y_plate**2) < (y_plate.max()/nFWHMy*nFWHMy_fit)**2)[0]
    
    phaseMap = phase_diff_2d[index_y[0]:index_y[-1]+1, index_x[0]:index_x[-1]+1]
    phaseMap -= phaseMap.min()
    
    # 2nd order polynomial fit
    ny, nx = phaseMap.shape
    fit_x = x_plate[index_x[0]:index_x[-1]+1]
    fit_y = y_plate[index_y[0]:index_y[-1]+1]
    soln = polyfit2d(fit_y, fit_x, phaseMap, kx=2, ky=2, order=2)
    
    xx, yy = np.meshgrid(x_plate, y_plate)
    fitted_surf = np.polynomial.polynomial.polyval2d(xx, yy, soln[0].reshape(3,3))
    
    delta, beta = load_CXRO_data(E0)
    plateThickness = (phase_diff_2d-fitted_surf) * 1239.84193e-9/E0 / (delta*2*np.pi)
    
    # plot
    plt.figure(figsize=(22,5))
    plt.subplot(1,4,1)
    plt.imshow(fitted_surf, cmap='jet',
               extent = [x_plate.min(), x_plate.max(), y_plate.max(), y_plate.min()])
    plt.colorbar(); plt.xlabel('2nd order fit'); plt.title('2nd order fit')
    plt.xlabel('x (mm)'); plt.ylabel('y (mm)'); plt.axis('tight')
    
    plt.subplot(1,4,2)
    plt.imshow(plateThickness*1e6, cmap='jet',
               extent = [x_plate.min(), x_plate.max(), y_plate.max(), y_plate.min()])
    plt.colorbar(); plt.xlabel('plate thickness (um)'); plt.title('plate thickness')
    plt.xlabel('x (mm)'); plt.ylabel('y (mm)'); plt.axis('tight')
    
    plt.subplot(1,4,3)
    plt.plot(x_plate, plateThickness[int(plateThickness.shape[0]/2)] * 1e6)
    plt.xlabel('x (mm)'); plt.ylabel('thickness (um)'); plt.title('lineout_x')
    
    plt.subplot(1,4,4)
    plt.plot(y_plate, plateThickness[:,int(plateThickness.shape[1]/2)] * 1e6)
    plt.xlabel('y (mm)'); plt.ylabel('thickness (um)'); plt.title('lineout_y')
    
    return plateThickness


''' calculates where a function crosses 0 in a given range '''
def find_zero(x, y, direction=None, x_i=0, x_f=1e5):
    # this function searches the solution of y(x) = 0 in range (x_i, x_f)
    # direction is 'increase' or 'decrease'
    y = y[np.intersect1d(np.where(x>=x_i), np.where(x<=x_f))]
    x = x[np.intersect1d(np.where(x>=x_i), np.where(x<=x_f))]

    # if linearly decrease:
    if np.median(y[:5]) >= np.median(y[-5:]):
        indices = np.where(y>=0)[0]
    else:
        indices = np.where(y<=0)[0]
    if direction == 'increase':
        indices = np.where(y<=0)[0]
    elif direction == 'decrease':
        indices = np.where(y>=0)[0]
    index = indices[-1]

    slope = (y[index+1] - y[index])/(x[index+1] - x[index])
    result = x[index] - y[index]/slope
    print('left {}, right {}, result {}'.format(x[index], x[index+1], result))
    return result