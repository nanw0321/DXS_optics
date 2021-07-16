class PhasePlate:
    """
    Class to represent arbitrary phase plate.

    Attributes
    ----------
    name: str
        Name of the device (e.g. CRL1)
    material: str
        Lens material. Currently only Be is implemented but may add CVD diamond in the future.
        Looks up downloaded data from CXRO.
    dx: float
        PhasePlate de-centering along beam's x-axis.
    dy: float
        PhasePlate de-centering along beam's y-axis.
    z: float
        z location of PhasePlate along beamline.
    energy: (N,) ndarray
        List of photon energies from CXRO file (eV).
    delta: (N,) ndarray
        Real part of index of refraction. n = 1 - delta + 1j * beta
    beta: (N,) ndarray
        Imaginary part of index of refraction. n = 1 - delta + 1j * beta
    """

    def __init__(self, name, plateThickness=None, x_mir=None, y_mir=None, E0=None, material='Be', z=0, dx=0, orientation=0):
        """
        Method to create a PhasePlate object.
        :param name: str
            Name of the device (e.g. CRL1)
        :param plateThickness: float
            Thickness of phase plate in [m]. Can be 1D (central line) or 2D.
        :param x_mir: float
            PhasePlate size in x.
        :param y_mir: float
            PhasePlate size in y.
        Looks up downloaded data from CXRO.
        :param z: float
            z location of lenses along beamline.
        :param dx: float
            Lens de-centering along beam's x-axis.
        :param dy: float
            Lens de-centering along beam's y-axis.
        """

        # set some attributes
        self.name = name
        self.plateThickness = plateThickness
        self.x_mir = x_mir
        self.y_mir = y_mir
        self.E0 = E0
        self.material = material
        self.dx = dx
        self.z = z
        self.global_x = 0
        self.global_y = 0
        self.orientation = orientation
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

        if self.orientation == 0:
            beamx = beam.x
            beamz = beam.zx
            beamc = beam.cx
            Ms = np.shape(self.plateThickness)[0]
            xs = np.linspace(-Ms/2, Ms/2 -1, Ms) * self.x_mir
        else:
            beamx = beam.y
            beamz = beam.zy
            beamc = beam.cy
            xs = self.y_mir
            Ns = np.shape(self.plateThickness)[1]
            xs = np.linspace(-Ns/2, Ns/2 -1, Ns) * self.y_mir

        # interpolate to find index of refraction at beam's energy
        delta = np.interp(beam.photonEnergy, self.energy, self.delta)
        beta = np.interp(beam.photonEnergy, self.energy, self.beta)

        # get shape of phase corrector input
        self.plate_shape = np.shape(self.plateThickness)

        # assume this is the central line of phase plate along the long axis if only 1D
        if np.size(self.plate_shape) == 1:
            central_line = self.plateThickness
        else:
            # get the central line from the 2D phase plate
            Ms = self.plate_shape[0]
            Ns = self.plate_shape[1]
            if self.orientation == 0:
                central_line = self.plateThickness[np.int(Ms/2)]
            else:
                central_line = self.plateThickness[:, np.int(Ns/2)]

        # 1D interpolation onto beam coordinates
        thickness = np.interp(beamx - self.dx/2, xs, central_line)

        # phase plate aperture
        mask = (((beamx - self.dx) ** 2) < (xs.max() / 2) ** 2).astype(float)

        phase = -beam.k0 * delta * thickness

        # transmission based on beta and thickness profile
        # phase shift at center of beam
        phase_shift = np.interp(beamc, beamx, thickness)*delta*2*np.pi/beam.lambda0

        transmission = np.exp(-beam.k0 * beta * thickness) * np.exp(1j * phase) * mask

        if self.orientation == 0:
            beam.wavex *= transmission
        else:
            beam.wavey *= transmission

    def propagate(self, beam):
        """
        Method to propagate beam through PhasePlate. Calls multiply.
        :param beam: Beam
            Beam object to propagate through PhasePlate. Beam is modified by this method.
        :return: None
        """
        self.multiply(beam)
