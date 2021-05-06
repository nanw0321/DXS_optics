import numpy as np
import os

bl = 'HRM'       # beamline to simulate: HHLM, HRM, HHLM_HRM
var_t = 400      # pulse duration [fs]
err_E = 0.       # error in photon energy [eV]
err_f0 = 0.      # error in collimation distance [m]
err_f1 = 0.      # error in crystal-lens distance [m]
C = 'C1'         # which crystal to add errors to
err_delta = 0.   # alignment error [rad]
err_miscut = 0.  # asymmetry angle error [rad]
if_log = 1       # linear or log scale for plotting
if_slice = 0     # x/y=0 slice or full for plotting
if_close = 1     # slit open or closed

# generate scripts
os.system(
    "python Gen_job_script.py --bl {} --var_t {} --err_E {} --err_f0 {} --err_f1 {} --C {} --err_delta {} --err_miscut {} --if_log {} --if_slice {} --if_close {}".format(
    bl, var_t, err_E, err_f0, err_f1, C, err_delta, err_miscut, if_log, if_slice, if_close))
