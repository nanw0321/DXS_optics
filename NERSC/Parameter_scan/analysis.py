import numpy as np
import os

var_t = 400.
bl_name = 'HRM'
err_name = '_f1'
err_f1_list = np.linspace(-1e-2, 1e-2, 50)

##### sampling parameters
sigT = var_t * 1e-15/2.355
d_slit = 7e-6

t_res =sigT*2.355/10              # time sampling resolution [s]
ev_window = 4/t_res *1e-15        # total energy window [eV]
ev_res = min(ev_window/800, 1e-3) # energy sampling resolution [eV]
nz = 2*round(ev_window/ev_res/2)
t_window = t_res*nz               # total time window [s]

# I/O
dir_output = 'output/'
dir_bl_out = dir_output+'{}/'.format(bl_name)    # HHLM, HRM, HHLM_HRM
dir_dur_out = dir_bl_out+'{}fs{}/'.format(round(var_t,1), err_name)    # pulse duration (100.0fs/)

diagnostic_names_open =  ['6_dur_out[fs]', '7_ptilt_x_out[fs_um]']
diagnostic_names_close = ['0_axis_ev_in[eV]', '1_int_ev_in[a.u.]', '2_axis_ev_out[eV]', '3_int_ev_out[a.u.]', '4_cent_E_out[eV]', '5_bw_out[eV]' ]

diagnostics_open = []
diagnostics_close = []

for job_num, err_f1 in enumerate(err_f1_list):
    job_num+=1
    err_val_name = '_{}m'.format(round(err_f1,3))
    dir_plot = dir_dur_out+'{}fs_{}meV_job{}{}/'.format(round(t_window*1e15,1), round(ev_window*1e3,1), job_num, err_val_name)    # sampling range + error value
    dir_close = dir_plot+'diagnostics_closed/'
    dir_open = dir_plot+'diagnostics_open/'

    for name in diagnostic_names_open:
        file_name = dir_open + name + '.npy'
        if not os.path.exists(file_name):
            print(file_name)
        else:
            with open(file_name, 'rb') as f:
                diagnostic = np.load(f)
                print(name, diagnostic.shape)
                diagnostics_open.append(diagnostic)

    for name in diagnostic_names_close:
        file_name = dir_close + name + '.npy'
        if not os.path.exists(file_name):
            print(file_name)
        else:
            with open(file_name, 'rb') as f:
                diagnostic = np.load(f)
                print(name, diagnostic.shape)
                diagnostics_close.append(diagnostic)

# for i, diagnostic in enumerate(diagnostics_open):
#     print(diagnostic_names_open[i], ':', diagnostic.shape)

# for i, diagnostic in enumerate(diagnostics_close):
#     print(diagnostic_names_close[i], ':', diagnostic.shape)