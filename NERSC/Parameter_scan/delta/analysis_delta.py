import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')   # allows plot without X11 forwarding

var_t = 100.
bl_name = 'HRM'
err_name = '_delta_C2'
err_delta_list = np.linspace(-1000, 1000, 50)

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

x = []
for job_num, err_delta in enumerate(err_delta_list):
    job_num+=1
    if job_num == 16: continue
    err_val_name = '_{}nrad'.format(round(err_delta,3))
    dir_plot = dir_dur_out+'{}fs_{}meV_job{}{}/'.format(round(t_window*1e15,1), round(ev_window*1e3,1), job_num, err_val_name)    # sampling range + error value
    dir_close = dir_plot+'diagnostics_closed/'
    dir_open = dir_plot+'diagnostics_open/'
    err_1 = 0; err_2 = 0
    for name in diagnostic_names_open:
        file_name = dir_open + name + '.out'
        if not os.path.exists(file_name):
            print(file_name)
            err_1 += 1
        else:
            diagnostic_open = np.loadtxt(file_name)
            diagnostics_open.append(diagnostic_open)

    for name in diagnostic_names_close:
        file_name = dir_close + name + '.out'
        if not os.path.exists(file_name):
            print(file_name)
            err_2 += 1
        else:
            diagnostic_close = np.loadtxt(file_name)
            diagnostics_close.append(diagnostic_close)

    if err_1+err_2 == 0:
        x.append(err_delta_list[job_num-1]/1e3)
    # if err_1 == 0 and err_2 != 1:
    #     diagnostics_close = diagnostics_close[:-len(diagnostic_names_close)]
    # if err_1 !=0 and err_2 == 0:
    #     diagnostics_open = diagnostics_open[:-len(diagnostic_names_open)]
    # if err_1 !=0 and err_2 != 0:
    #     diagnostics_close = diagnostics_close[:-len(diagnostic_names_close)]
    #     diagnostics_open = diagnostics_open[:-len(diagnostic_names_open)]

# plot
x = np.array(x) # [mrad]
xlabel = 'error delta {} (mrad)'.format(err_name[-2:])
fig, axes = plt.subplots(nrows=1, ncols=5, sharex=True, figsize=(27,5))
# plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

# 1. pulse duration
duration = np.array(diagnostics_open[0::len(diagnostic_names_open)])
ax1 = axes[0]
ax1.plot(x, duration)
ax1.set_xlabel(xlabel)
ax1.set_ylabel('fs')
ax1.set_title('pulse duration')

# 2. pulse front tilt
ptilt = np.array(diagnostics_open[1::len(diagnostic_names_open)])
ax2 = axes[1]
ax2.plot(x, ptilt)
ax2.set_xlabel(xlabel)
ax2.set_ylabel('fs/um')
ax2.set_title('pulse-front tilt')

# 3. bandwidth
bw = np.array(diagnostics_close[5::len(diagnostic_names_close)])*1e3
ax3 = axes[2]
ax3.plot(x, bw)
ax3.set_xlabel(xlabel)
ax3.set_ylabel('meV')
ax3.set_title('bandwidth')

# 4. througput
axis_ev_in = np.array(diagnostics_close[0::len(diagnostic_names_close)])    # 50x800
axis_ev_out = np.array(diagnostics_close[2::len(diagnostic_names_close)])
int_ev_in = np.array(diagnostics_close[1::len(diagnostic_names_close)])
int_ev_out = np.array(diagnostics_close[3::len(diagnostic_names_close)])

axis_ev = axis_ev_in[0]
throughput = int_ev_out/int_ev_in
cent_E = (axis_ev[0]+axis_ev[-1])/2
axis_ev = (axis_ev - cent_E)*1e3

ax4 = axes[3]
im = ax4.pcolormesh(x, axis_ev, throughput.T, cmap='jet')
plt.colorbar(im, ax=ax4)
ax4.set_xlabel(xlabel)
ax4.set_ylabel('meV around {}eV'.format(cent_E,2))
ax4.set_title('spectral response')

# 5. central energy
Ec = np.array(diagnostics_close[4::len(diagnostic_names_close)])
ax5 = axes[4]
ax5.plot(x, (Ec-cent_E)*1e3)
ax5.set_xlabel(xlabel)
ax5.set_ylabel('meV around {}eV'.format(round(cent_E,2)))
ax5.set_title('central energy')

plt.savefig('scan{}.png'.format(err_name))