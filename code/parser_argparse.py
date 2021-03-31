##### parse parameter
# python test.py -var_t 100 -err_E 10 -err_f0 1e-1 -err_f1 1e-1 -C C1 -err_delta 1e-6 -err_miscut 1e-6
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--var_t', type=float, required=False,
        help='variable: pulse_duration [fs]')
    parser.add_argument('--err_E', type=float, required=False,
        help='error: photon energy [meV]')
    parser.add_argument('--err_f0', type=float, required=False,
        help='error: collimation lens focal distance [m]')
    parser.add_argument('--err_f1', type=float, required=False,
        help='error: crystal-lens distance [m]')
    parser.add_argument('--C', type=str, required=False,
        help='crystal name [string]')
    parser.add_argument('--err_delta', type=float, required=False,
        help='error: crystal alignment [urad]')
    parser.add_argument('--err_miscut', type=float, required=False,
        help='error: crystal miscut [urad]')
    parser.add_argument('--if_log', type=int, required=False,
        help='condition: log scale plots')
    parser.add_argument('--if_slice', type=int, required=False,
        help='condition: y=0 slice plots')
    parser.add_argument('--if_close', type=int, required=False,
        help='condition: closed slit')
    

    args = parser.parse_args()

    var_t = args.var_t
    err_E = args.err_E
    err_f0 = args.err_f0
    err_f1 = args.err_f1

    cname = args.C
    err_delta = args.err_delta
    err_miscut = args.err_miscut
    
    if_log = args.if_log
    if_slice = args.if_slice
    if_close = args.if_close

    err_name = ''
    err_val_name = ''

    if var_t is None: var_t = 400.

    if err_E is None: err_E = 0.
    else: err_name += '_E'; err_val_name += '_{}meV'.format(round(err_E,3))

    if err_f0 is None: err_f0 = 0.
    else: err_name += '_f0'; err_val_name += '_{}m'.format(round(err_f0,3))

    if err_f1 is None: err_f1 = 0.
    else: err_name += '_f1'; err_val_name += '_{}m'.format(round(err_f1,3))
        
    if cname is None: cname = 'HHLM1'
        
    if err_delta is None: err_delta = 0.
    else: err_name += '_delta_{}'.format(cname); err_val_name += '_{}urad'.format(round(err_delta,3))
        
    if err_miscut is None: err_miscut = 0.
    else: err_name += '_miscut_{}'.format(cname); err_val_name += '_{}urad'.format(round(err_miscut,3))

    if if_log is None: if_log = 1
    if if_slice is None: if_slice = 1
    if if_close is None: if_close = 0

    print('arguments: {}'.format([var_t, err_E, err_f0, err_f1, err_delta, err_miscut]))
    print('error type: {}, error value: {}'.format(err_name, err_val_name))
    print('if_log {}, if_close {}, if_slice {}'.format(if_log, if_close, if_slice))