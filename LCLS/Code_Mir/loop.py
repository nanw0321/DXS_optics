from Functions import *
import scipy.interpolate as interpolate

def load_data(num):
    dir_profile = 'DXS_Xtals_FEA4WFS/'
    coords = np.loadtxt(dir_profile+'Nlist_%d.txt' %num, skiprows=12)
    data = np.loadtxt(dir_profile+'Uy_list_%d.txt' % num, skiprows=0)

    x = coords[:,1]
    y = coords[:,2]
    z = coords[:,3]
    dy = data[:,1]
    print(data.shape)
    # plt.figure()
    # plt.plot(x)
    # plt.figure()
    # plt.plot(z)

    xx = np.linspace(np.min(x), np.max(x), 1024)
    zz = np.linspace(np.min(z), np.max(z), 1024)
    xx,zz = np.meshgrid(xx,zz)

    # f = interpolate.interp2d(z,x,dy)
    dy2 = interpolate.griddata((x,z),dy,(xx,zz), fill_value=0)

    dy_symmetrize = np.concatenate((dy2,np.flipud(dy2)),axis=0)
    dy_symmetrize = np.concatenate((dy_symmetrize,np.fliplr(dy_symmetrize)),axis=1)

    xx2 = np.linspace(np.min(x), -np.min(x), 2048)
    zz2 = np.linspace(np.min(z), -np.min(z), 2048)
    xx2,zz2 = np.meshgrid(xx2,zz2)
    
    return dy_symmetrize, xx2, zz2

def plot_data(num, dy, xx, zz):
    plt.figure(figsize=(7 ,6))
    ax_im = plt.subplot2grid((6,7), (1,0), colspan=5, rowspan=3)
    ax_x = plt.subplot2grid((6,7), (4,0), colspan=5, rowspan=2)
    ax_y = plt.subplot2grid((6,7), (1,5), colspan=2, rowspan=3)
    ax_colorbar = plt.subplot2grid((6,7), (4,5), colspan=1, rowspan=2)
    im_profile = ax_im.imshow(dy*1e6+125, extent=(-7.5,7.5,-.75,.75),aspect=(20/4))
    ax_im.set_xlabel('Tangential coordinates (mm)')
    ax_im.set_ylabel('Sagittal coordinates (mm)')
    ax_x.plot(xx[0,:],dy[1024,:]*1e6+125)
    ax_x.set_ylabel('Height (nm)')
    ax_x.grid()
    ax_y.plot(dy[:,1024]*1e6+125,zz[:,0])
    ax_y.set_xlabel('Height (nm)')
    ax_y.grid()
    plt.colorbar(im_profile, cax=ax_colorbar, label='Height error (nm)')
    ax_im.set_title('Case %d' % num)

#### Load crystal shapeError
crystals = np.array(['HHLM1', 'HHLM2', 'HHLM3', 'HHLM4', 'C1', 'C2'])    # list of crystals
n_crys = len(crystals)

option = None    # SASE, reduced SASE, self-seeded
if_error = np.array([1, 0, 0, 0, 0, 0])    # if add shapeError to each crystal
height_factor = 1.0
slit_factor = 1.0

x_FWHMs = [14456.77e-6, 14459.84e-6, 14452.53e-6, 1312.47e-6, 1317.63e-6, 321.91e-6]    # beam size at the output of each crystal
y_FWHMs = [706.86e-6, 707.26e-6, 707.94e-6, 708.32e-6, 718.2e-6, 719.41e-6]    # beam size at the output of each crystal
beta0s = [0.6441930561063302, 0.7492421541043375, 0.7492421541043375, 0.05102436691260805, 0.7492421541043375, 0.23448464050278972]    # output angle
projected_x_FWHMs = np.divide(x_FWHMs, np.sin(beta0s))    # beam FWHMs projected on crystal surface

for i_iter in range(8):
    ''' choose between SASE, SASE w/ reduced power for HRM, or self-seeded '''
    if option is not None:
        options = [option for i in range(n_crys)]    # SASE, SASE w/ reduced power for HRM, self-seeded
        if_error = [1 for i in range(n_crys)]    # turn on all shapeError
    else:
        options = [1, 3, 3, 3, 3, 3]    # cases for largest error
        
    # some other special cases
    if i_iter == 0:
        # simulate when shapeError is reduced by half, HHLM2 only
        option = None
        options[1] = 3
        if_error = [0, 1, 0, 0, 0, 0]
        height_factor = 0.5

    if i_iter == 1:
        # simulate if shapeError at 16W can be corrected for HHLM2
        option = None
        options[1] = 1
        if_error = [0, 1, 0, 0, 0, 0]
        height_factor = 1.0

    if i_iter >1:
        # simulate if our simulation is sensitive to different slit width
        option = None
        options = [1, 3, 3, 3, 3, 3]
        if_error = [0 for i in range(n_crys)]; if_error[i_iter-2] = 1
        slit_factor = 2.0

        print('shapeError on/off',if_error)

    crystal_name = crystals[np.where(np.array(if_error)==1)[0]]
    crystal_option = np.array(options)[np.where(np.array(if_error)==1)[0]]
    print(options, crystal_name, crystal_option)

    shapeErrors = [None for i in range(n_crys)]
    x_mirs      = [None for i in range(n_crys)]
    y_mirs      = [None for i in range(n_crys)]
    l_crystal   = [10e-2 for i in range(n_crys)]    # HRM crystal lengths (long)
    w_crystal   = [1e-2 for i in range(n_crys)]    # HRM crystal widths (short)

    for i, crystal in enumerate(crystals):
        if if_error[i] == 1:
            num = (options[i]-1)*n_crys + i+1
            dy, xx, zz = load_data(num)
            shapeErrors[i] = dy * 1e6
            x_mirs[i] = xx[0]
            y_mirs[i] = zz[:,0]

    ##### Padding
    plt.figure(figsize=(int(5*n_crys)+2, 10))
    for i in range(n_crys):
        if if_error[i] == 1:
            plt.subplot(2,n_crys, int(i+1))
            plt.plot(x_mirs[i], shapeErrors[i][1024,:], label=crystals[i])
            plt.legend()
            plt.xlabel('x (mm)')
            plt.ylabel('y (nm)')
            plt.title('{}, {}mm FWHM'.format(crystals[i], round(projected_x_FWHMs[i]*1e3,2)))

            shape = shapeErrors[i]

            shapeError_x = (x_mirs[i].max() - x_mirs[i].min())*1e-3; dx = (x_mirs[i][1]-x_mirs[i][0])*1e-3
            shapeError_y = (y_mirs[i].max() - y_mirs[i].min())*1e-3; dy = (y_mirs[i][1]-y_mirs[i][0])*1e-3

            l_crystal[i] = shapeError_x
            w_crystal[i] = shapeError_y

            size_x_padded = np.int(l_crystal[i]/dx) + 1; size_x_padded = size_x_padded + size_x_padded%2
            size_y_padded = np.int(w_crystal[i]/dy) + 1; size_y_padded = size_y_padded + size_y_padded%2

            shapeError_new = np.ones((size_x_padded, size_y_padded)) * shape.max()
            roi_x_i = int((size_x_padded - shape.shape[0])/2)
            roi_x_f = roi_x_i + shape.shape[0]
            roi_y_i = int((size_y_padded - shape.shape[1])/2)
            roi_y_f = roi_y_i + shape.shape[1]

            shapeError_new[roi_x_i:roi_x_f, roi_y_i:roi_y_f] = shape.T
            shapeError_new = shapeError_new.T * height_factor
            
            plt.subplot(2,n_crys, n_crys+(int(i+1)))
            plt.imshow(shapeError_new, cmap='jet')
            plt.colorbar()
            plt.title('padded')
            plt.axis('tight')

            shapeErrors[i] = shapeError_new

    ''' Beam parameters '''
    N = 1024            # number of sampling points
    E0 = 9481            # photon energy [eV]
    tau = 400            # pulse duration [fs]
    window = 50000        # total simulation time window [fs]

    # parameter dictionary. z_source is in LCLS coordinates (20 meters upstream of undulator exit)
    beam_params = {
        'photonEnergy': E0,
        'N': N,
        'sigma_x': 23e-6,
        'sigma_y': 23e-6,
        'rangeFactor': 5,
        'scaleFactor': 10,
        'z_source': 650-26
    }

    ##### I/O
    dir_out = 'output/'; make_dir(dir_out)
    dir_ev = dir_out + '{}eV_{}x/'.format(E0, round(slit_factor,1)); make_dir(dir_ev)

    if option is not None:
        dir_plot = dir_ev + '{}/'.format(option)
        if height_factor != 1:
            dir_plot = dir_ev + '{}_{}/'.format(option, height_factor)
        make_dir(dir_plot)
    else:
        dir_plot = dir_ev + '{}/'.format(crystal_name)
        if height_factor != 1:
            dir_plot = dir_ev +'{}_{}/'.format(crystal_name, height_factor)
        make_dir(dir_plot)

    print('\nCWD:',dir_plot)

    def define_beamline(
        E0, slit_width, m2_p=175.332, shapeErrors=[None for i in range(6)]):
        
        z_s = 650

        # viewing point - upstream of monochromator
        im_input = optics.PPM('im_input', z=184+z_s, FOV=5e-3, N=256)
        
        ## Telescope
        m1 = optics.CurvedMirror('M1',p=185,q=-58,length=1,z=185+z_s,alpha=2.65e-3)
        im_after_T1 = optics.PPM('im_after_T1', z=m1.z+.01, FOV=5e-3, N=256)
        
        m2 = optics.CurvedMirror('M2',p=m2_p,q=1e5,length=1,z=300+z_s,alpha=2.65e-3,orientation=2)
        im_after_T2 = optics.PPM('im_after_T2', z=m2.z+.01, FOV=5e-3, N=256)
        
        ## HHLM
        asym = np.deg2rad(17)
        hhlm1 = optics.Crystal('HHLM1', hkl=[2,2,0], length=l_crystal[0], width=w_crystal[0], z=305+z_s, alphaAsym=-asym, E0=E0, 
                               orientation=0, pol='s', shapeError=shapeErrors[0])
        im_after_HHLM1 = optics.PPM('im_after_HHLM1', FOV=30e-3,N=256,z=hhlm1.z+.01)
        
        hhlm2 = optics.Crystal('HHLM2', hkl=[4,4,0], length=l_crystal[1], width=w_crystal[1], z=hhlm1.z+139e-3, E0=E0,
                               orientation=2, pol='s', shapeError=shapeErrors[1])
        im_after_HHLM2 = optics.PPM('im_after_HHLM2', FOV=30e-3,N=256,z=hhlm2.z+.01)
        
        hhlm3 = optics.Crystal('HHLM3', hkl=[4,4,0], length=l_crystal[2], width=w_crystal[2], z=hhlm1.z+361e-3, E0=E0,
                               orientation=0, pol='s', shapeError=shapeErrors[2])
        im_after_HHLM3 = optics.PPM('im_after_HHLM3', FOV=30e-3,N=256,z=hhlm3.z+.01)
        
        hhlm4 = optics.Crystal('HHLM4', hkl=[2,2,0], length=l_crystal[3], width=w_crystal[3], z=hhlm1.z+.5, alphaAsym=asym, E0=E0,
                               orientation=2, pol='s', shapeError=shapeErrors[3])
        im_after_HHLM4 = optics.PPM('im_after_HHLM4', FOV=30e-3,N=256,z=hhlm4.z+.01)
        
        ## HRM
        alphaAsym = np.deg2rad(29.5)    # crystal asymmetry angle
        hkl = [4,4,0]    # crystal plane
        f1 = 10.    # crystal-lens/mirror distance
        f2 = 10.    # lens/mirror focal distance
        
        # first crystal: symmetric reflection
        crystal1 = optics.Crystal('C1', hkl=hkl, length=l_crystal[4], width=w_crystal[4], z=z_s+310, E0=E0, alphaAsym=0, 
                                  orientation=0, pol='s', delta=0.e-6, shapeError=shapeErrors[4])
        im_after_C1 = optics.PPM('im_after_C1', z=crystal1.z+.01, FOV=5e-3, N=256)

        # second crystal: asymmetric reflection, orientation flipped relative to crystal1
        crystal2 = optics.Crystal('C2', hkl=hkl, length=l_crystal[5], width=w_crystal[5], z=crystal1.z+.05, E0=E0,alphaAsym=alphaAsym, 
                                  orientation=2, pol='s', delta=0e-6, shapeError=shapeErrors[5])

        im_after_C2 = optics.PPM('im_after_C2', z=crystal2.z+.01, FOV=5e-3, N=256)
        
        # printing crystal incidence and reflection angles for confirmation
        print('crystal 2 incidence angle: {:.2f} degrees'.format(crystal2.alpha*180/np.pi))
        print('crystal 2 exit angle: {:.2f} degrees'.format(crystal2.beta0*180/np.pi))

        # elliptical mirror with ~10 meter focal length
        im_before_MIR1 = optics.PPM('im_before_MIR1', z=crystal2.z + f1 - .1, FOV=2e-3, N=256)
        mir1 = optics.CurvedMirror('mir1', z=crystal2.z+f1, p=1e5, q=f2, length=1.0, width=5e-3, alpha=3.6e-3, orientation=0)
        im_after_MIR1 = optics.PPM('im_after_MIR1', z=mir1.z+.1, FOV=5e-3, N=256)

        # slit at focus
        slit = optics.Slit('Slit', z=mir1.z+f2, x_width=slit_width, y_width=2e-3)
        print('slit width: {} um'.format(slit.x_width*1e6))

        # viewing point at focus
        im_focus = optics.PPM('im_focus', z=mir1.z+f2 + 1e-3, FOV=1e-3, N=1024)

        # elliptical mirror with ~10 meter focal length, for collimation
        im_before_MIR2 = optics.PPM('im_before_MIR2', z=mir1.z+2*f2 -.1, FOV=2e-3, N=256)
        mir2 = optics.CurvedMirror('mir2', z=mir1.z+2*f2, p=f2, q=1e5, length=1.0, width=5e-3, alpha=3.6e-3, orientation=2)
        im_after_MIR2 = optics.PPM('im_after_MIR2', z=mir2.z+.1, FOV=2e-3, N=256)
        
        # third crystal, symmetric reflection, same orientation as crystal2
    #     crystal3 = optics.Crystal('c3', hkl=hkl, length=10e-2, width=10e-3, z=mir2.z+d23, E0=E0,alphaAsym=0, orientation=2,
    #                              asym_type='emergence',pol='s')
        crystal3 = optics.Crystal('C3', hkl=hkl, length=10e-2, width=3e-3, z=mir2.z+f1, E0=E0,alphaAsym=-alphaAsym,
                                  orientation=2, pol='s')
        im_after_C3 = optics.PPM('im_after_C3', z=crystal3.z+.01, FOV=5e-3, N=256)

        # fourth crystal, asymmetric reflection, same orientation as crystal1
    #     crystal4 = optics.Crystal('c4', hkl=hkl, length=10e-2, width=10e-3, z=mir2.z+d23 + (f1-d23)*np.cos(crystal1.beta0*2), E0=E0,alphaAsym=-alphaAsym, 
    #                               asym_type='emergence', orientation=0,pol='s')
        crystal4 = optics.Crystal('C4', hkl=hkl, length=10e-2, width=3e-3, z=crystal3.z+.05, E0=E0,alphaAsym=0, 
                                  orientation=0, pol='s')
        im_out = optics.PPM('im_out', z=crystal4.z+.1, FOV=5e-3, N=256)

        # list of devices to propagate through
        telescope_devices = [im_input, m1, im_after_T1, m2, im_after_T2]
        HHLM_devices = [hhlm1, im_after_HHLM1, hhlm2, im_after_HHLM2, hhlm3, im_after_HHLM3, hhlm4, im_after_HHLM4]
        HRM_devices = [crystal1,im_after_C1, crystal2,im_after_C2, im_before_MIR1,mir1,im_after_MIR1, slit,im_focus,
                   im_before_MIR2,mir2,im_after_MIR2, crystal3,im_after_C3, crystal4,im_out]
        
        all_devices = telescope_devices + HHLM_devices + HRM_devices
        mono_beamline = beamline.Beamline(all_devices, ordered=True)
        return mono_beamline


    #### calculate slit width
    import sys, os
    sys.stdout = open(os.devnull, 'w')

    def calc_slit_width(inbeam, lmbd, foc):
        # inbeam is the beamsize (2x FWHM) after Crystal 2
        # lmbd is wavelength of photon
        # asym is asymmetry angle of Crystal 2
        # foc is Mir1/Lens1 focal distance
        # returns mono focus size in [m] at the Fourier plane
        slt = 2*np.log(2)*lmbd*foc/np.pi/inbeam
        return 2.11 * slt    # 2x FWHM with some error (~factor of 5%)

    ##### define beamline
    if_close = 1    # if close slit

    if if_close == 1:
        # propagate to calculate beam size after crystal 2
        beam_params['photonEnergy'] = E0
        mono_beamline = define_beamline(E0,200e-6)
        b1 = beam.Beam(beam_params=beam_params)
        # pulse = beam.Pulse(beam_params=beam_params, tau=tau, time_window=window)
        pulse = beam.Pulse(beam_params=beam_params, unit_spectrum=True, spectral_width=40e-3, N=200)
        screen_names = ['im_after_C2']
        pulse.propagate(beamline=mono_beamline, screen_names=screen_names)
        _,_, wx, wy, _,_ = get_width(pulse, 'im_after_C2')
        slit_width = calc_slit_width(wx*2, b1.lambda0, mono_beamline.mir1.q)*slit_factor
    else:
        slit_width = 200e-6

    mono_beamline = define_beamline(E0,slit_width)

    ##### prepare for loops
    n_loop = 31
    spectral_width = 40e-3
    N = 200
    # m2_p = np.linspace(-15,15,n_loop) + 175.332
    m2_p = np.linspace(150,400, n_loop)

    pulse_duration0 = np.zeros(n_loop); pulse_duration1 = np.zeros(n_loop)
    pulse_tilt0 = np.zeros(n_loop); pulse_tilt1 = np.zeros(n_loop)
    spatial_chirp0 = np.zeros(n_loop); spatial_chirp1 = np.zeros(n_loop)
    energy_bandwidth0 = np.zeros(n_loop); energy_bandwidth1 = np.zeros(n_loop)

    fname_perfect = dir_plot+'diagnostics_perfect.h5'
    fname_error = dir_plot+'diagnostics_error.h5'
    if if_close == 1:
        fname_perfect = dir_plot+'diagnostics_perfect_close.h5'
        fname_error = dir_plot+'diagnostics_error_close.h5'

    ##### loop, Perfect crystals
    for i in range(n_loop):
        # reset central photon energy
        beam_params['photonEnergy'] = E0

        # initialize pulse. tau sets transform-limited pulse duration in fs (FWHM). Time window is in fs.
        pulse = beam.Pulse(beam_params=beam_params, unit_spectrum=True, spectral_width=spectral_width, N=N)
        mono_beamline = define_beamline(E0,slit_width, m2_p=m2_p[i],shapeErrors=[None for i in range(n_crys)])
        
        pulse.propagate(beamline=mono_beamline,screen_names=['im_out'])
        c, pulse_duration0[i] = pulse.pulse_duration('im_out')
        pulse_tilt0[i] = pulse.pulsefront_tilt('im_out')
        spatial_chirp0[i] = pulse.spatial_chirp('im_out')
        energy_bandwidth0[i] = pulse.pulse_bandwidth('im_out')
        
    with h5py.File(fname_perfect, 'w') as f:
        f.create_dataset('m2_p', data=m2_p)
        f.create_dataset('pulse_duration', data=pulse_duration0)
        f.create_dataset('tilt', data=pulse_tilt0)
        f.create_dataset('chirp', data=spatial_chirp0)
        f.create_dataset('bandwidth', data=energy_bandwidth0)

    ##### loop, With Error
    for i in range(n_loop):
        # reset central photon energy
        beam_params['photonEnergy'] = E0

        # initialize pulse. tau sets transform-limited pulse duration in fs (FWHM). Time window is in fs.
        pulse = beam.Pulse(beam_params=beam_params, unit_spectrum=True, spectral_width=spectral_width, N=N)
        mono_beamline = define_beamline(E0,slit_width, m2_p=m2_p[i],shapeErrors=shapeErrors)
        
        pulse.propagate(beamline=mono_beamline,screen_names=['im_out'])
        c, pulse_duration1[i] = pulse.pulse_duration('im_out')
        pulse_tilt1[i] = pulse.pulsefront_tilt('im_out')
        spatial_chirp1[i] = pulse.spatial_chirp('im_out')
        energy_bandwidth1[i] = pulse.pulse_bandwidth('im_out')

    with h5py.File(fname_error, 'w') as f:
        f.create_dataset('m2_p', data=m2_p)
        f.create_dataset('pulse_duration', data=pulse_duration1)
        f.create_dataset('tilt', data=pulse_tilt1)
        f.create_dataset('chirp', data=spatial_chirp1)
        f.create_dataset('bandwidth', data=energy_bandwidth1)

    sys.stdout = sys.__stdout__
    
    #### plot
    fname_perfect = dir_plot+'diagnostics_perfect.h5'
    fname_error = dir_plot+'diagnostics_error.h5'
    if if_close == 1:
        fname_perfect = dir_plot+'diagnostics_perfect_close.h5'
        fname_error = dir_plot+'diagnostics_error_close.h5'

    with h5py.File(fname_perfect, 'r') as f:
        m2_p = f['m2_p'][:]
        pulse_duration0 = f['pulse_duration'][:]
        pulse_tilt0 = f['tilt'][:]
        spatial_chirp0 = f['chirp'][:]
        energy_bandwidth0 = f['bandwidth'][:]
        
    with h5py.File(fname_error, 'r') as f:
        m2_p = f['m2_p'][:]
        pulse_duration1 = f['pulse_duration'][:]
        pulse_tilt1 = f['tilt'][:]
        spatial_chirp1 = f['chirp'][:]
        energy_bandwidth1 = f['bandwidth'][:]

    if option == 1: suptitle = 'Full shapeError, SASE'
    elif option == 2: suptitle = 'Full shapeError, Reduced SASE for HRM'
    elif option == 3: suptitle = 'Full shapeError, Self-seeded'
    else: suptitle = '{} only shapeError, option {}'.format(crystal_name, crystal_option)

    suptitle = suptitle+', M2 p {}m'.format(175.332)
    if if_close == 1:
        suptitle = suptitle + ', slit_width {}um'.format(round(slit_width*1e6,3))

    plt.figure(figsize=(22,5))
    plt.suptitle(suptitle)
    plt.subplot(1,4,1)
    plt.plot(m2_p, pulse_duration0,label='perfect')
    plt.plot(m2_p, pulse_duration1,label='heat bump')
    plt.xlabel('M2 p (m)')
    plt.ylabel('pulse duration (fs)')
    # plt.ylim([0, 5000])
    plt.grid()
    plt.legend()
    plt.tight_layout()

    plt.subplot(1,4,2)
    plt.plot(m2_p, pulse_tilt0,label='perfect')
    plt.plot(m2_p, pulse_tilt1,label='heat bump')
    plt.xlabel('M2 p (m)')
    plt.ylabel(r'pulsefront tilt (fs/$\mu$m)')
    plt.grid()
    plt.legend()
    plt.tight_layout()

    plt.subplot(1,4,3)
    plt.plot(m2_p, spatial_chirp0*1000,label='perfect')
    plt.plot(m2_p, spatial_chirp1*1000,label='heat bump')
    plt.xlabel('M2 p (m)')
    plt.ylabel(r'spatial chirp (meV/$\mu$m)')
    plt.grid()
    plt.legend()
    plt.tight_layout()

    plt.subplot(1,4,4)
    plt.plot(m2_p, energy_bandwidth0*1000,label='perfect')
    plt.plot(m2_p, energy_bandwidth1*1000,label='heat bump')
    plt.xlabel('M2 p (m)')
    plt.ylabel('bandwidth (meV)')
    plt.grid()
    plt.legend()
    plt.tight_layout()

    if option is not None:
        pltname = 'option_{}_correction'.format(option)
    else:
        pltname = '{}_correction_option_{}'.format(crystal_name, crystal_option)

    if if_close == 1:
        pltname += '_close'
    plt.savefig(dir_ev+'{}.png'.format(pltname))
    print(dir_plot)

