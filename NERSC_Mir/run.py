import numpy as np
import os

bl = 'HRM'       # beamline to simulate: HHLM, HRM, HHLM_HRM
var_t = 20      # pulse duration [fs]
err_E = 0.       # error in photon energy [eV]
err_f0 = 0.      # error in collimation distance [m]
err_f1 = 0.      # error in crystal-lens distance [m]
C = 'C1'         # which crystal to add errors to
err_delta = 0.   # alignment error [rad]
err_miscut = 0.  # asymmetry angle error [rad]
if_log = 0       # linear or log scale for plotting
if_slice = 0     # x/y=0 slice or full for plotting
if_close = 0     # slit open or closed

# generate scripts
for job_num, var_t in enumerate([20, 40, 100, 200, 400]):
    os.system(
        "python Gen_job_script.py --bl {} --var_t {} --err_E {} --err_f0 {} --err_f1 {} --C {} --err_delta {} --err_miscut {} --if_log {} --if_slice {} --if_close {} --job {}".format(
            bl, var_t, err_E, err_f0, err_f1, C, err_delta, err_miscut, if_log, if_slice, 0, job_num))
    os.system(
        "python Gen_job_script.py --bl {} --var_t {} --err_E {} --err_f0 {} --err_f1 {} --C {} --err_delta {} --err_miscut {} --if_log {} --if_slice {} --if_close {} --job {}".format(
            bl, var_t, err_E, err_f0, err_f1, C, err_delta, err_miscut, if_log, if_slice, 1, job_num))

