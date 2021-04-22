import numpy as np
import os

bl = 'HRM'
var_t = 400
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
err_E_list = np.linspace(-1000, 1000, 2)	# [meV]
for job_num, err_E in enumerate(err_E_list):
    if job_num >0: continue
    os.system(
        "python Gen_job_script.py --bl {} --var_t {} --err_E {} --err_f0 {} --err_f1 {} --C {} --err_delta {} --err_miscut {} --if_log {} --if_slice {} --if_close {} --job {}".format(
            bl, var_t, err_E, err_f0, err_f1, C, err_delta, err_miscut, if_log, if_slice, 0, job_num+1))
    os.system(
        "python Gen_job_script.py --bl {} --var_t {} --err_E {} --err_f0 {} --err_f1 {} --C {} --err_delta {} --err_miscut {} --if_log {} --if_slice {} --if_close {} --job {}".format(
            bl, var_t, err_E, err_f0, err_f1, C, err_delta, err_miscut, if_log, if_slice, 1, job_num+1))