


''' get info '''
def print_oe_pos(oe):
    print('{}, x:{}, y:{}, z:{}'.format(oe.name, oe.global_x, oe.global_y, oe.z))
    return oe.global_x, oe.global_y, oe.z

def get_width(pulse, image_name):
    # minima and maxima of the field of view (in microns) for imshow extent
    minx = np.round(np.min(pulse.x[image_name]) * 1e6)
    maxx = np.round(np.max(pulse.x[image_name]) * 1e6)
    miny = np.round(np.min(pulse.y[image_name]) * 1e6)
    maxy = np.round(np.max(pulse.y[image_name]) * 1e6)

    # calculate the profile
    profile = np.sum(np.abs(pulse.energy_stacks[image_name])**2, axis=2)
    x_lineout = np.sum(profile, axis=0)
    y_lineout = np.sum(profile, axis=1)

    cx, cy, wx, wy, fwx_guess, fwy_guess = pulse.beam_analysis(pulse.x[image_name],pulse.y[image_name],
                                                                       x_lineout, y_lineout)

    return cx, cy, wx, wy, fwx_guess, fwy_guess

def get_pulse(pulse, image_name, x_pos=0, y_pos=0, shift=None):
    minx = np.round(np.min(pulse.x[image_name]) * 1e6)
    maxx = np.round(np.max(pulse.x[image_name]) * 1e6)
    miny = np.round(np.min(pulse.y[image_name]) * 1e6)
    maxy = np.round(np.max(pulse.y[image_name]) * 1e6)

    # get number of pixels
    M = pulse.x[image_name].size
    N = pulse.y[image_name].size

    # calculate pixel sizes (microns)
    dx = (maxx - minx) / M
    dy = (maxy - miny) / N

    # calculate indices for the desired location
    x_index = int((x_pos - minx) / dx)
    y_index = int((y_pos - miny) / dy)

    # calculate temporal intensity
    y_data = np.abs(pulse.time_stacks[image_name][y_index, x_index, :]) ** 2

    shift = -pulse.t_axis[np.argmax(y_data)]

    # coarse shift for fitting
    if shift is not None:
        y_data = np.roll(y_data, int(shift/pulse.deltaT))

    # get gaussian stats
    centroid, sx = Util.gaussian_stats(pulse.t_axis, y_data)
    fwhm = int(sx * 2.355)

    # gaussian fit
    gauss_plot = Util.fit_gaussian(pulse.t_axis, centroid, sx)

    # shift again using fit result
    shift = -centroid
    if shift is not None:
        y_data = np.roll(y_data, int(shift/pulse.deltaT))
        gauss_plot = np.roll(gauss_plot, int(shift/pulse.deltaT))
        
    # [fs], normalized intensity [simulated], [Gaussian Fit]
    return pulse.t_axis, y_data/np.max(y_data), gauss_plot

def get_spectrum(pulse, image_name, x_pos=0, y_pos=0, integrated=False):
    minx = np.round(np.min(pulse.x[image_name]) * 1e6)
    maxx = np.round(np.max(pulse.x[image_name]) * 1e6)
    miny = np.round(np.min(pulse.y[image_name]) * 1e6)
    maxy = np.round(np.max(pulse.y[image_name]) * 1e6)

    # get number of pixels
    M = pulse.x[image_name].size
    N = pulse.y[image_name].size

    # calculate pixel sizes (microns)
    dx = (maxx - minx) / M
    dy = (maxy - miny) / N

    # calculate indices for the desired location
    x_index = int((x_pos - minx) / dx)
    y_index = int((y_pos - miny) / dy)

    # calculate spectral intensity
    if integrated:
        y_data = np.sum(np.abs(pulse.energy_stacks[image_name])**2, axis=(0,1))
    else:
        y_data = np.abs(pulse.energy_stacks[image_name][y_index,x_index,:])**2

    # get gaussian stats
    centroid, sx = Util.gaussian_stats(pulse.energy, y_data)
    fwhm = sx * 2.355

    # gaussian fit to plot
    gauss_plot = Util.fit_gaussian(pulse.energy, centroid, sx)

    # change label depending on bandwidth
    if fwhm >= 1:
        width_label = '%.1f eV FWHM' % fwhm
    elif fwhm > 1e-3:
        width_label = '%.1f meV FHWM' % (fwhm * 1e3)
    else:
        width_label = u'%.1f \u03BCeV FWHM' % (fwhm * 1e6)
    
    # [eV], normalized intensity [simulated], [Gaussian Fit]
    return pulse.energy - pulse.E0, y_data/np.max(y_data), gauss_plot


''' plot '''
def plot_phase_comparison(ppm_p, ppm_e, ppm_c, x_range, y_range, image_name=None):
    # crop roi based on beam width
    index_x = np.where((ppm_p.x **2) < (x_range ** 2))[0]
    index_y = np.where((ppm_p.y **2) < (y_range ** 2))[0]
    axis_x = np.linspace(-x_range, x_range, index_x[-1] - index_x[0] + 1)*1e3
    axis_y = np.linspace(-y_range, y_range, index_y[-1] - index_y[0] + 1)*1e3

    # 1D
    phase_p_x = ppm_p.x_phase[index_x[0]:index_x[-1]+1]
    phase_p_y = ppm_p.y_phase[index_y[0]:index_y[-1]+1]
    phase_e_x = ppm_e.x_phase[index_x[0]:index_x[-1]+1]
    phase_e_y = ppm_e.y_phase[index_y[0]:index_y[-1]+1]
    phase_c_x = ppm_c.x_phase[index_x[0]:index_x[-1]+1]
    phase_c_y = ppm_c.y_phase[index_y[0]:index_y[-1]+1]

    plt.figure(figsize=(30,20))
    plt.subplot(2,2,1)
    plt.plot(axis_x, phase_p_x - np.median(phase_p_x), label='w/o')
    plt.plot(axis_x, phase_e_x - np.median(phase_e_x), label='w/')
    plt.plot(axis_x, phase_c_x - np.median(phase_c_x), label='corrected')
    plt.title('phase x, {}'.format(image_name)); plt.xlabel('x (mm)')
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(axis_y, phase_p_y - np.median(phase_p_y), label='w/o')
    plt.plot(axis_y, phase_e_y - np.median(phase_e_y), label='w/')
    plt.plot(axis_y, phase_c_y - np.median(phase_c_y), label='corrected')
    plt.title('phase y'); plt.xlabel('y (mm)')
    plt.legend()
