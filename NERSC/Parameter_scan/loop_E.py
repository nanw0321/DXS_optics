import numpy as np
import os

bl = 'HRM'
var_t = 100
err_E = 0.
err_f0 = 0.
err_f1 = 0.
C = 'C1'
err_delta = 0.
err_miscut = 0.
if_log = 1
if_slice = 0
if_close = 0

# generate scripts and submit job
err_E_list = np.linspace(-5, 5, 5)	# [meV]
for job_num, err_E in enumerate(err_E_list):
    os.system(
        "python Gen_job_script.py --bl {} --var_t {} --err_E {} --err_f0 {} --err_f1 {} --C {} --err_delta {} --err_miscut {} --if_log {} --if_slice {} --if_close {} --job {}".format(
            bl, var_t, err_E, err_f0, err_f1, C, err_delta, err_miscut, if_log, if_slice, if_close, job_num+1))
