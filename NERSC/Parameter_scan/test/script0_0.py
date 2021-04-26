from Module_diagnostic_functions import *
from Module_propagation_functions import *

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
    
    return wfr
    # t0 = time(); print('Resizing in frequency domain: ', end='')
    # srwlpy.ResizeElecField(wfr, 'f', [0, 1., z_scaling]); print('done in', round(time() - t0, 3), 's')

    # # diagnostics - input
    # diagnostics, diagnostics_names = diagnose_input(wfr)

    # # CC1
    # label = 'before C1'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    # bl = set_optics_CRL0_before_C1(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    # plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    # label = 'before C2'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    # bl = set_optics_C1_before_C2(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    # plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    # label = 'before CRL1'; print('Propagating to {}: '.format(label), end=''); i_plot+=1; t0 = time()
    # bl = set_optics_C2_before_CRL1(v); srwlpy.PropagElecField(wfr, bl); print('done in', round(time() - t0, 3), 's\n')
    # plot_wfr_diagnostic(wfr, label=label, dir_plot=dir_plot, i=i_plot, if_log=if_log, if_slice=if_slice, i_start=i_start)

    # # diagnostics - output
    # diagnostics, diagnostics_names = diagnose_output(wfr, diagnostics=diagnostics, diagnostics_names=diagnostics_names)

    # # I/O
    # if if_close == 1:
    #     dir_diagnostics = dir_plot+'diagnostics_closed/'; mkdir(dir_diagnostics)
    # else:
    #     dir_diagnostics = dir_plot+'diagnostics_open/'; mkdir(dir_diagnostics)
    # save_diagnostics(dir_diagnostics, diagnostics, diagnostics_names)

    # print('\nbeam diagnostics written to {}'.format(dir_diagnostics))

if __name__ == '__main__':
    main_HRM('varParam0.pkl', z_scaling=6.0, if_close=0, dir_plot='output/HRM/400.0fs/4000.0fs_100.0meV_job0/', if_log=1, if_slice=0, i_start=10)