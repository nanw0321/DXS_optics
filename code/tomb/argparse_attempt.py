##### diagnostic
from Diagnostic_functions import *
import argparse

parser = argparse.ArgumentParser(description='list of sampling parameters')
parser.add_argument("-X", type = float, help = 'xRange resize factor')
parser.add_argument("-x", type = float, help = 'xRes resize factor')
parser.add_argument("-Y", type = float, help = 'yRange resize factor')
parser.add_argument("-y", type = float, help = 'yRes resize factor')
parser.add_argument("-s", type = float, help = 'slit width (um)')
parser.add_argument("-f", type= float, help = 'CRL focal length (m)')
parser.add_argument("-t", type= float, help = 'pulse duration (fs)')
args = parser.parse_args()

if args.X    is not None: xRange = args.X
else: xRange = 20
if args.x    is not None: xRes = args.x
else: xRes = 4
if args.Y   is not None: yRange = args.Y
else: yRange = 1
if args.y is not None: yRes = args.y
else: yRes = 2
if args.s    is not None: d_slit = args.s/1e6
else: d_slit = 1e-3
if args.f    is not None: fCRL1 = args.f
else: fCRL1 = 10.
if args.t    is not None: sigT = args.t/1e15
else: sigT = 400e-15

print(__name__)
print(xRange, xRes, yRange, yRes, d_slit, fCRL1, sigT)