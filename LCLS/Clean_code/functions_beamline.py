import time, h5py, os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from lcls_beamline_toolbox.xraybeamline2d import beam1d as beam, optics1d as optics, beamline1d as beamline
from lcls_beamline_toolbox.xraybeamline2d.util import Util

''' define beamline 1: telescope mirrors T1, T2 '''
def define_Telescope(E0, z_s=650.0, m_alpha=2.65e-3, m2_z=300.0,
                     m1_p=185.0, m1_q=-58.0, m2_p=None, m2_q=None,
                     aperture_size=1.0):
    '''
    defines the Telescope optics
    E0: photon energy [eV]
    m_alpha: grazing angle [rad]
    m1_p: mirror 1 source distance [m]
    m1_q: mirror 1 image distance [m], negative value indicates diverging mirror
    m2_p: mirror 2 source distance [m]
    m2_q: mirror 2 image distance [m]
    aperture size: # +- nxFWHM aperture size for telescope mirrors [m]

    returns optics
    '''
    
    ## Telescope
    m1 = optics.CurvedMirror('M1', p=m1_p, q=m1_q, length=1*aperture_size, z=m1_p+z_s, alpha=m_alpha, orientation=0)
    if m2_p is None:
        m2_p = m2_z - m1_p - m1_q
    if m2_q is None:
        m2_q = 1e5    # fix the image distance to be 100km so the output is almost parallel
    m2 = optics.CurvedMirror('M2', p=m2_p, q=m2_q, length=1*aperture_size, z=m2_z+z_s, alpha=m_alpha, orientation=2)
    
    im_after_T1 = optics.PPM('im_after_T1', z=m1.z+.01, FOV=5e-3, N=512)
    im_after_T2 = optics.PPM('im_after_T2', z=m2.z+.01, FOV=5e-3, N=512)

    Telescope_devices = [m1, im_after_T1, m2, im_after_T2]

    return Telescope_devices


''' define beamline 2: HHLM (high heatload mono) '''
def define_HHLM_2DCM(E0, z_s=650.0,
                     HHLM_offset=20e-3,
                     pair_distance=200e-3,
                     hkl1 = [1,1,1], alphaAsym1 = 9.0,
                     hkl2 = [1,1,1], alphaAsym2 = 0.0,
                     shapeErrors=[None for i in range(6)],
                     l_crystal=[1e-1 for i in range(6)],
                     w_crystal = [5e-3 for i in range(6)]):
    
    '''
    defines the HHLM optics for the 2DCM setup (1-1-2-2)
    E0: photon energy [eV]
    HHLM_offset: beam offset between crystal 1 and 2
    pair_distance: distance between crytal 2 and 3
    hkl: crystal reflection surface indices for pair 1 and pair 2
    alphaAsym: asymmetry angle [degree], negative value means beam footprint increases after reflection
    shapeErrors: crytal shapeError as loaded from Lin's profiles

    returns optics
    '''

    ## HHLM
    asym1 = np.deg2rad(alphaAsym1)
    asym2 = np.deg2rad(alphaAsym2)
    hhlm1 = optics.Crystal('HHLM1', hkl=hkl1, length=l_crystal[0], width=w_crystal[0],
                           z=305+z_s, alphaAsym=-asym1, E0=E0, orientation=0, pol='s',
                           shapeError=shapeErrors[0])
    d12 = HHLM_offset/np.tan(2*hhlm1.bragg)
    hhlm2 = optics.Crystal('HHLM2', hkl=hkl1, length=l_crystal[1], width=w_crystal[1],
                           z=hhlm1.z+d12, alphaAsym=asym1, E0=E0, orientation=2, pol='s',
                           shapeError=shapeErrors[1])

    hhlm3 = optics.Crystal('HHLM3', hkl=hkl2, length=l_crystal[2], width=w_crystal[2],
                           z=hhlm2.z+pair_distance, alphaAsym=-asym2, E0=E0, orientation=2, pol='s',
                           shapeError=shapeErrors[2])
    d34 = HHLM_offset/np.tan(2*hhlm3.bragg)
    hhlm4 = optics.Crystal('HHLM4', hkl=hkl2, length=l_crystal[3], width=w_crystal[3],
                           z=hhlm3.z+d34, alphaAsym=asym2, E0=E0, orientation=0, pol='s',
                           shapeError=shapeErrors[3])

    im_after_HHLM1 = optics.PPM('im_after_HHLM1', FOV=2.5e-2,N=512,z=hhlm1.z+d12/100)
    im_after_HHLM2 = optics.PPM('im_after_HHLM2', FOV=5e-3,N=512,z=hhlm2.z+1e-3)
    im_after_HHLM3 = optics.PPM('im_after_HHLM3', FOV=2.5e-2,N=512,z=hhlm3.z+d34/100)
    im_after_HHLM4 = optics.PPM('im_after_HHLM4', FOV=5e-3,N=512,z=hhlm4.z+1e-3)

    HHLM_devices = [hhlm1, im_after_HHLM1, hhlm2, im_after_HHLM2, hhlm3, im_after_HHLM3, hhlm4, im_after_HHLM4]

    return HHLM_devices

def define_HHLM_Zigzag(E0, z_s=650.0,
                     HHLM_offset=20e-3,
                     pair_distance=200e-3,
                     hkl1 = [1,1,1], alphaAsym1 = 9.0,
                     hkl2 = [1,1,1], alphaAsym2 = 0.0,
                     shapeErrors=[None for i in range(6)],
                     l_crystal=[1e-1 for i in range(6)],
                     w_crystal = [5e-3 for i in range(6)]):
    '''
    defines the HHLM optics for the zigzag setup (1-2-2-1)
    E0: photon energy [eV]
    HHLM_offset: beam offset between crystal 1 and 2
    pair_distance: distance between crytal 2 and 3
    hkl: crystal reflection surface indices for pair 1 and pair 2
    alphaAsym: asymmetry angle [degree], negative value means beam footprint increases after reflection
    shapeErrors: crytal shapeError as loaded from Lin's profiles

    returns optics
    '''

    ## HHLM
    asym1 = np.deg2rad(alphaAsym1)
    asym2 = np.deg2rad(alphaAsym2)
    hhlm1 = optics.Crystal('HHLM1', hkl=hkl1, length=l_crystal[0], width=w_crystal[0],
                           z=305+z_s, alphaAsym=-asym1, E0=E0, orientation=0, pol='s',
                           shapeError=shapeErrors[0])
    d12 = HHLM_offset/np.tan(2*hhlm1.bragg)
    hhlm2 = optics.Crystal('HHLM2', hkl=hkl2, length=l_crystal[1], width=w_crystal[1],
                           z=hhlm1.z+d12, alphaAsym=-asym2, E0=E0,orientation=2, pol='s',
                           shapeError=shapeErrors[1])
    if hhlm2.bragg-hhlm1.bragg < np.pi/4: 
        d23 = pair_distance
    else: 
        d23 = -pair_distance
    hhlm3 = optics.Crystal('HHLM3', hkl=hkl2, length=l_crystal[2], width=w_crystal[2],
                           z=hhlm2.z+d23, alphaAsym=asym2, E0=E0,orientation=0, pol='s',
                           shapeError=shapeErrors[2])
    d34 = np.abs(pair_distance*np.tan(2*(hhlm2.bragg-hhlm1.bragg))-HHLM_offset)/np.tan(2*hhlm1.bragg)
    hhlm4 = optics.Crystal('HHLM4', hkl=hkl1, length=l_crystal[3], width=w_crystal[3],
                           z=hhlm3.z+d34, alphaAsym=asym1, E0=E0,orientation=2, pol='s',
                           shapeError=shapeErrors[3])
    
    im_after_HHLM1 = optics.PPM('im_after_HHLM1', FOV=2.5e-2,N=512,z=hhlm1.z+d12/100)
    im_after_HHLM2 = optics.PPM('im_after_HHLM2', FOV=2.5e-2,N=512,z=hhlm2.z+np.sign(d23)*1e-3)
    im_after_HHLM3 = optics.PPM('im_after_HHLM3', FOV=5e-3,N=512,z=hhlm3.z+d34/100)
    im_after_HHLM4 = optics.PPM('im_after_HHLM4', FOV=5e-3,N=512,z=hhlm4.z+1e-3)

    HHLM_devices = [hhlm1, im_after_HHLM1, hhlm2, im_after_HHLM2, hhlm3, im_after_HHLM3, hhlm4, im_after_HHLM4]

    return HHLM_devices


''' define beamline 3: HRM (high res mono) '''
def define_HRM(E0, z_s=650.0,
               f1=10., f2=10., slit_width=3e-6,
               hkl=[4,4,0], alphaAsym=15.0,
               shapeErrors=[None for i in range(6)],
               l_crystal=[1e-1 for i in range(6)],
               w_crystal = [5e-3 for i in range(6)]):
    '''
    defines the HRM optics
    E0: photon energy [eV]
    f1: crystal-lens/mirror distance [m]
    f2: lens/mirror focal distance [m]
    slit_width: width of slit for monochromatization [m]
    hkl: crystal reflection surface indices for pair 1 and pair 2
    alphaAsym: asymmetry angle [degree], negative value means beam footprint increases after reflection
    shapeErrors: crytal shapeError as loaded from Lin's profiles

    returns optics
    '''

    ## HRM
    asym3 = np.deg2rad(alphaAsym)
    crystal1 = optics.Crystal('C1', hkl=hkl, length=l_crystal[4], width=w_crystal[4],
                              z=z_s+310, E0=E0, alphaAsym=0, orientation=0, pol='s',
                              shapeError=shapeErrors[4])
    d12 = 0.02/np.tan(2*crystal1.bragg)
    crystal2 = optics.Crystal('C2', hkl=hkl, length=l_crystal[5], width=w_crystal[5],
                              z=crystal1.z+d12, E0=E0,alphaAsym=asym3, orientation=2, pol='s',
                              shapeError=shapeErrors[5])
    
    mir1 = optics.CurvedMirror('mir1', z=crystal2.z+f1, p=1e5, q=f2, length=1.0, width=5e-3, alpha=3.6e-3, orientation=0)
    
    # slit at focus
    slit = optics.Slit('Slit', z=mir1.z+f2*np.cos(2*mir1.alpha), x_width=slit_width, y_width=2e-3)
    print('slit width: {} um'.format(slit.x_width*1e6))
    
    mir2 = optics.CurvedMirror('mir2', z=mir1.z+2*f2*np.cos(2*mir1.alpha), p=f2, q=1e5, length=1.0, width=5e-3, alpha=3.6e-3, orientation=2)
    
    crystal3 = optics.Crystal('C3', hkl=hkl, length=1e-1, width=5e-3,
                              z=mir2.z+f1, E0=E0,alphaAsym=-asym3, orientation=2, pol='s')
    d34 = (0.02 + 2*f2*np.sin(2*mir1.alpha))/np.tan(2*crystal3.bragg)

    crystal4 = optics.Crystal('C4', hkl=hkl, length=1e-1, width=5e-3,
                              z=crystal3.z+d34, E0=E0,alphaAsym=0, orientation=0, pol='s')
    
    im_after_C1    = optics.PPM('im_after_C1',    z=crystal1.z+d12/100, FOV=5e-3, N=512)
    im_after_C2    = optics.PPM('im_after_C2',    z=crystal2.z+1e-3, FOV=5e-3, N=512)
    im_before_MIR1 = optics.PPM('im_before_MIR1', z=mir1.z-1e-3,     FOV=5e-3, N=512)
    im_after_MIR1  = optics.PPM('im_after_MIR1',  z=mir1.z+1e-3,     FOV=5e-3, N=512)
    im_focus       = optics.PPM('im_focus',       z=slit.z+1e-3,     FOV=50e-6, N=512)
    im_before_MIR2 = optics.PPM('im_before_MIR2', z=mir2.z-1e-3,     FOV=5e-3, N=512)
    im_after_MIR2  = optics.PPM('im_after_MIR2',  z=mir2.z+1e-3,     FOV=5e-3, N=512)
    im_after_C3    = optics.PPM('im_after_C3',    z=crystal3.z+d34/100, FOV=5e-3, N=512)
    im_after_C4    = optics.PPM('im_after_C4',    z=crystal4.z+1e-3, FOV=5e-3, N=512)
    
    HRM_devices = [crystal1,im_after_C1, crystal2,im_after_C2, im_before_MIR1,mir1,im_after_MIR1, slit,im_focus,
                   im_before_MIR2,mir2,im_after_MIR2, crystal3,im_after_C3, crystal4,im_after_C4]
    
    return HRM_devices


''' define overall beamline '''
def define_beamline_normal(E0, z_s=650, m_alpha=2.65e-3, m2_z=300.0,
                           m1_p=None, m1_q=None, m2_p=None, m2_q=None,
                           aperture_size=1.0,

                           HHLM_type='2DCM',
                           HHLM_offset=20e-3,
                           pair_distance=200e-3,
                           hkl1 = [1,1,1], alphaAsym1 = 0.0,
                           hkl2 = [2,2,0], alphaAsym2 = 0.0,
                           l_crystal=[1e-1 for i in range(6)],
                           w_crystal = [5e-3 for i in range(6)],
                           
                           f1=10.0, f2=10.0, slit_width=3e-6,
                           hkl3 = [5,5,5], alphaAsym3 = 15.0,
                           shapeErrors=[None for i in range(6)],
                           
                           FOV1 = 5e-3, N1 = 512,
                           FOV2 = 5e-3, N2 = 8192,
                           plate_position = 'HHLM'):
    
    # partial beamlines    
    Telescope_devices = define_Telescope(E0, z_s=z_s, m_alpha=m_alpha, m2_z=m2_z,
                                         m1_p=m1_p, m1_q=m1_q, m2_p=m2_p, m2_q=m2_q,
                                         aperture_size=aperture_size)
    
    if HHLM_type == '2DCM':
        HHLM_devices = define_HHLM_2DCM(E0, z_s=z_s,
                                        HHLM_offset=HHLM_offset,
                                        pair_distance=pair_distance,
                                        hkl1=hkl1, alphaAsym1=alphaAsym1,
                                        hkl2=hkl2, alphaAsym2=alphaAsym2,
                                        shapeErrors=shapeErrors,
                                        l_crystal=l_crystal,
                                        w_crystal=w_crystal)
    elif HHLM_type == 'Zigzag':
        HHLM_devices = define_HHLM_Zigzag(E0, z_s=z_s,
                                        HHLM_offset=HHLM_offset,
                                        pair_distance=pair_distance,
                                        hkl1=hkl1, alphaAsym1=alphaAsym1,
                                        hkl2=hkl2, alphaAsym2=alphaAsym2,
                                        shapeErrors=shapeErrors,
                                        l_crystal=l_crystal,
                                        w_crystal=w_crystal)
    
    HRM_devices = define_HRM(E0, z_s=z_s,
                             f1=f1, f2=f2, slit_width=slit_width,
                             hkl=hkl3, alphaAsym=alphaAsym3,
                             shapeErrors=shapeErrors,
                             l_crystal=l_crystal,
                             w_crystal=w_crystal)
    
    # viewing points
    im_input = optics.PPM('im_input', z=184+z_s, FOV=FOV1, N=N1)
    im_output = optics.PPM('im_output', FOV=FOV1,N=N1,z=HRM_devices[-1].z+1.1e-3)
    
    # putting it together
    if plate_position == 'HHLM':
        im_plate = optics.PPM('im_plate', FOV=FOV2, N=N2, z=HHLM_devices[-1].z + 1e-3)
        all_devices = [im_input] + Telescope_devices + HHLM_devices + [im_plate] + HRM_devices + [im_output]
    elif plate_position == 'HRM':
        im_plate = optics.PPM('im_plate', FOV=FOV2, N=N1, z=HHLM_devices[-1].z + 1e-3)
        im_plate2 = optics.PPM('im_plate2', FOV=FOV2, N=N2, z=HRM_devices[-1].z + 1e-3)
        all_devices = [im_input] + Telescope_devices + HHLM_devices + [im_plate] + HRM_devices + [im_plate2, im_output]
    
    mono_beamline = beamline.Beamline(all_devices, ordered=True)
    return all_devices, mono_beamline


''' phase plate '''
class PhasePlate:
    """
    Attributes
    ----------
    name: str
        Name of the device (e.g. CRL1)
    plateThickness: float
        Thickness profile of the phase plate. (meters)
    x_plate: float
        Phase plate size in x. (meters)
    y_plate: float
        Phase plate size in y. (meters)
    E0: float or None
        photon energy in eV for calculating the corresponding phase difference of a given thickness
    material: str
        Phase plate material. Currently only Be is implemented but may add CVD diamond in the future.
        Looks up downloaded data from CXRO.
    dx: float
        Phase plate de-centering along beam's x-axis.
    dy: float
        Phase plate de-centering along beam's y-axis.
    z: float
        z location of phase plate along beamline.
    energy: (N,) ndarray
        List of photon energies from CXRO file (eV).
    delta: (N,) ndarray
        Real part of index of refraction. n = 1 - delta + 1j * beta
    beta: (N,) ndarray
        Imaginary part of index of refraction. n = 1 - delta + 1j * beta
    """
    
    def __init__(self, name, plateThickness=None, x_plate=None, y_plate=None, E0=None, material='Be', z=0, dx=0, dy=0):
        """
        Method to create a PhasePlate object.
        :param name: str
            Name of the device (e.g. Phase1)
        :param plateThickness: float
            Thickness profile of the phase plate. (meters)
        :x_plate: float
            Phase plate size in x. (meters)
        :y_plate: float
            Phase plate size in y. (meters)
        :param E0: float
            photon energy for calculating radius of curvature for a given focal length (eV)
        :param material: str
            Lens material. Currently only Be is implemented but may add CVD diamond in the future.
            Looks up downloaded data from CXRO.
        :param z: float
            z location of lenses along beamline.
        :param dx, dy: float
            PhasePlate de-centering along beam's x,y-axis.
        :param orientation: int
            Whether or not this is a horizontal or vertical lens (0 for horizontal, 1 for vertical).
        """
        
        # set some attributes
        self.name = name
        self.plateThickness = plateThickness
        self.x_plate = x_plate
        self.y_plate = y_plate
        self.E0 = E0
        self.material = material
        self.dx = dx
        self.dy = dy
        self.z = z
        self.global_x = 0
        self.global_y = 0
        self.azimuth = 0
        self.elevation = 0
        self.xhat = None
        self.yhat = None
        self.zhat = None

        # get file name of CXRO data
        filename = os.path.join('../../LCLS/lcls_beamline_toolbox-beta/lcls_beamline_toolbox/xraybeamline2d/cxro_data/Be.csv')

        # load in CXRO data
        cxro_data = np.genfromtxt(filename, delimiter=',')
        self.energy = cxro_data[:, 0]
        self.delta = cxro_data[:, 1]
        self.beta = cxro_data[:, 2]

    def multiply(self, beam):
        """
        Method to propagate beam through PhasePlate
        :param beam: Beam
            Beam object to propagate through PhasePlate. Beam is modified by this method.
        :return: None
        """
        
        # get shape of phase plate thickness
        self.plate_shape = np.shape(self.plateThickness)
        Ms = self.plate_shape[0]; Ns = self.plate_shape[1]
        
        
        beamx = beam.x
        beamy = beam.y
        central_line_x = self.plateThickness[np.int(Ms/2)]
        central_line_y = self.plateThickness[:, np.int(Ns/2)]
        
        Ms = central_line_x.size; Ns = central_line_y.size
        xs = np.linspace(-Ms/2, Ms/2 -1, Ms) * self.x_plate/Ms    # phase plate x coordinate
        ys = np.linspace(-Ns/2, Ns/2 -1, Ns) * self.y_plate/Ns    # phase plate y coordinate
        
        # interpolation onto beam coordinates
        thickness_x = np.interp(beamx - self.dx, xs, central_line_x, left=0, right=0)
        thickness_y = np.interp(beamy - self.dy, ys, central_line_y, left=0, right=0)

        # interpolate to find index of refraction at beam's energy
        delta = np.interp(beam.photonEnergy, self.energy, self.delta)
        beta = np.interp(beam.photonEnergy, self.energy, self.beta)
        phase_x = -beam.k0 * delta * thickness_x
        phase_y = -beam.k0 * delta * thickness_y
        
        # transmission based on beta and thickness profile
        mask_x = (((beamx - self.dx) ** 2) < (self.x_plate / 2) ** 2).astype(float)
        mask_y = (((beamy - self.dy) ** 2) < (self.y_plate / 2) ** 2).astype(float)
        transmission_x = np.exp(-beam.k0 * beta * thickness_x) * np.exp(1j * phase_x) * mask_x
        transmission_y = np.exp(-beam.k0 * beta * thickness_y) * np.exp(1j * phase_y) * mask_y
        
        beam.wavex *= transmission_x
        beam.wavey *= transmission_y

    def propagate(self, beam):
        """
        Method to propagate beam through PhasePlate. Calls multiply.
        :param beam: Beam
            Beam object to propagate through PhasePlate. Beam is modified by this method.
        :return: None
        """
        self.multiply(beam)


def define_beamline_phase(E0, z_s=650, m_alpha=2.65e-3, m2_z=300.0,
                           m1_p=None, m1_q=None, m2_p=None, m2_q=None,
                           aperture_size=1.0,

                           HHLM_type='2DCM',
                           HHLM_offset=20e-3,
                           pair_distance=200e-3,
                           hkl1 = [1,1,1], alphaAsym1 = 0.0,
                           hkl2 = [2,2,0], alphaAsym2 = 0.0,
                           l_crystal=[1e-1 for i in range(6)],
                           w_crystal = [5e-3 for i in range(6)],
                           
                           f1=10.0, f2=10.0, slit_width=3e-6,
                           hkl3 = [5,5,5], alphaAsym3 = 15.0,
                           shapeErrors=[None for i in range(6)],
                           
                           FOV1 = 5e-3, N1 = 512,
                           FOV2 = 5e-3, N2 = 8192,
                           plate_position = 'HHLM',
                           plate_length=1., plate_width=1.,
                           plateThickness = None):
    
    # partial beamlines    
    Telescope_devices = define_Telescope(E0, z_s=z_s, m_alpha=m_alpha, m2_z=m2_z,
                                         m1_p=m1_p, m1_q=m1_q, m2_p=m2_p, m2_q=m2_q,
                                         aperture_size=aperture_size)
    
    if HHLM_type == '2DCM':
        HHLM_devices = define_HHLM_2DCM(E0, z_s=z_s,
                                        HHLM_offset=HHLM_offset,
                                        pair_distance=pair_distance,
                                        hkl1=hkl1, alphaAsym1=alphaAsym1,
                                        hkl2=hkl2, alphaAsym2=alphaAsym2,
                                        shapeErrors=shapeErrors,
                                        l_crystal=l_crystal,
                                        w_crystal=w_crystal)
    elif HHLM_type == 'Zigzag':
        HHLM_devices = define_HHLM_Zigzag(E0, z_s=z_s,
                                        HHLM_offset=HHLM_offset,
                                        pair_distance=pair_distance,
                                        hkl1=hkl1, alphaAsym1=alphaAsym1,
                                        hkl2=hkl2, alphaAsym2=alphaAsym2,
                                        shapeErrors=shapeErrors,
                                        l_crystal=l_crystal,
                                        w_crystal=w_crystal)

    HRM_devices = define_HRM(E0, z_s=z_s,
                             f1=f1, f2=f2, slit_width=slit_width,
                             hkl=hkl3, alphaAsym=alphaAsym3,
                             shapeErrors=shapeErrors,
                             l_crystal=l_crystal,
                             w_crystal=w_crystal)
    
    # viewing points
    im_input = optics.PPM('im_input', z=184+z_s, FOV=FOV1, N=N1)
    im_output = optics.PPM('im_output', FOV=FOV1,N=N1,z=HRM_devices[-1].z+1.1e-3)
    
    # putting it together
    if plate_position == 'HHLM':
        Phase1 = PhasePlate('Phase1', plateThickness=plateThickness,
                            x_plate=plate_length, y_plate=plate_width,
                            E0=E0, z=HHLM_devices[-1].z + 1e-3)
        im_plate = optics.PPM('im_plate', FOV=FOV2, N=N2, z=HHLM_devices[-1].z + 1.1e-3)
        all_devices = [im_input] + Telescope_devices + HHLM_devices + [Phase1, im_plate] + HRM_devices + [im_output]
    elif plate_position == 'HRM':
        Phase1 = PhasePlate('Phase1', plateThickness=plateThickness,
                            x_plate=plate_length, y_plate=plate_width,
                            E0=E0, z=HHLM_devices[-1].z + 1e-3)
        im_plate = optics.PPM('im_plate', FOV=FOV2, N=N1, z=HHLM_devices[-1].z + 1e-3)
        im_plate2 = optics.PPM('im_plate2', FOV=FOV2, N=N2, z=HRM_devices[-1].z + 1.1e-3)
        all_devices = [im_input] + Telescope_devices + HHLM_devices + [im_plate] + HRM_devices + [Phase1, im_plate2, im_output]
    
    mono_beamline = beamline.Beamline(all_devices, ordered=True)
    return all_devices, mono_beamline


''' load crystal shapeErrors from file '''
def load_crystal_data(dir_profile, crystal_name, option,
                      nFWHMx=3, nFWHMx_profile=3,
                      nFWHMy=1, nFWHMy_profile=1):
    coords = np.loadtxt(dir_profile+'Nlist_{}_{}.txt'.format(option, crystal_name), skiprows=12)
    data = np.loadtxt(dir_profile+'Uy_list_{}_{}.txt'.format(option, crystal_name), skiprows=0)

    x = coords[:,1]
    y = coords[:,2]
    z = coords[:,3]
    dy = data[:,1]
    print(dy.shape)
    
    # intepolate onto the targeted FWHM range
    xx = np.linspace(np.min(x), np.max(x), 1024) / nFWHMx_profile * nFWHMx
    zz = np.linspace(np.min(z), np.max(z), 1024) / nFWHMy_profile * nFWHMy
    xx, zz = np.meshgrid(xx,zz)

    dy2 = interpolate.griddata((x,z), dy, (xx,zz))
    dy_symmetrize = np.concatenate((dy2, np.flipud(dy2)),axis=0)
    dy_symmetrize = np.concatenate((dy_symmetrize, np.fliplr(dy_symmetrize)),axis=1)
    
    xx2 = np.linspace(np.min(xx), -np.min(xx), 2048)
    zz2 = np.linspace(np.min(zz), -np.min(zz), 2048)
    xx2,zz2 = np.meshgrid(xx2,zz2)
    return dy_symmetrize, xx2, zz2


''' beamline optics adjustment '''
def crystal_alignment_error(devices, delta, crystal_name):
    # change crystal alignment angle delta
    for i, device in enumerate(devices):
        if device.name == crystal_name:
            devices[i].delta = delta

def crystal_miscut_error(devices, eta_err, crystal_name):
    # add crystal miscut error
    for i, device in enumerate(devices):
        if device.name == crystal_name:
            devices[i] = optics.Crystal(device.name, hkl=device.hkl, length=device.length, width=device.width,
                                        z=device.z, E0=device.E0, alphaAsym=device.alphaAsym+eta_err,
                                        orientation=device.orientation, pol=device.pol, delta=device.delta)

def oe_positional_shift(devices, shift, oe_name):
    # add positional shift along lab z axis
    for i, device in enumerate(devices):
        if device.name == oe_name:
            devices[i].z += shift

def lens_energyError(devices, E):
    # shift lens central energy
    for i, device in enumerate(devices):
        if device.name[:3] == 'crl':
            devices[i] = optics.CRL(device.name, z=device.z, E0=E, f=device.f, diameter=device.diameter)
