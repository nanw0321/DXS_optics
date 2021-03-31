import sys

err_name = ''
err_val_name = ''

try: var_t = float(sys.argv[1])
except: var_t = 400.

try: err_E = float(sys.argv[2])
except: err_E = 0.
if err_E != 0.: err_name += '_E'; err_val_name += '_{}meV'.format(round(err_E,3))

try: err_f0 = float(sys.argv[3])
except: err_f0 = 0.
if err_f0 != 0.: err_name += '_f0'; err_val_name += '_{}m'.format(round(err_f0,3))

try: err_f1 = float(sys.argv[4])
except: err_f1 = 0.
if err_f1 != 0.: err_name += '_f1'; err_val_name += '_{}m'.format(round(err_f1,3))

try: cname = str(sys.argv[5])
except: cname = 'HHLM1'

try: err_delta = float(sys.argv[6])
except: err_delta = 0.
if err_delta != 0.: err_name += '_delta_{}'.format(cname); err_val_name += '_{}urad'.format(round(err_delta,3))

try: err_miscut = float(sys.argv[7])
except: err_miscut = 0.
if err_miscut != 0.: err_name += '_miscut_{}'.format(cname); err_val_name += '_{}urad'.format(round(err_miscut,3))

try: if_log = int(sys.argv[8])
except: if_log = 1

try: if_close = int(sys.argv[9])
except: if_close = 0
    
try: if_slice = int(sys.argv[10])
except: if_slice = 1

print('arguments: {}'.format([var_t, err_E, err_f0, err_f1, cname, err_delta, err_miscut]))
print('error type: {}, error value: {}'.format(err_name, err_val_name))
print('if_log {}, if_close {}, if_slice {}'.format(if_log, if_close, if_slice))