##### parse parameter
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-var_t', type=float, required=True,
    help='variable: pulse_duration [fs]')
parser.add_argument('-var_E', type=float, required=True,
    help='variable: photon energy [eV]')
parser.add_argument('-err_f1', type=float, required=True,
    help='error: crystal-lens distance [m]')
parser.add_argument('-err_delta', type=float, required=True,
    help='error: crystal alignment [urad]')
parser.add_argument('-err_miscut', type=float, required=True,
    help='error: crystal miscut [urad]')
args = parser.parse_args()

var_t = args.var_t
var_E = args.var_E
err_f1 = args.err_f1
err_delta = args.err_delta
err_miscut = args.err_miscut

print('arguments: {}'.format([var_t, var_E, err_f1, err_delta, err_miscut]))


# Si 220
Si220 = srwlib.SRWLOptCryst(_d_sp=1.9201374688016222, _psi0r=-1.0873035035585694e-05, _psi0i=1.8438837339536554e-07,
                         _psi_hr=-6.610167427198717e-06, _psi_hi=1.7824173540780476e-07,
                         _psi_hbr=-6.610167427198717e-06, _psi_hbi=1.7824173540780476e-07,
                         _tc=0.01, _ang_as=0)
thetaB220 = Si220.get_ang_inc(_e=E)