from .simulation import start_simu, reg_wave, bretschneider_wave, decay_test
import COERbuoy.GUI as GUIServer
import argparse
from argparse import RawTextHelpFormatter as tf;

parser=argparse.ArgumentParser(description="The COERbuoy model",formatter_class=tf)
parser.add_argument("--regular_wave", nargs=4,type=str, metavar=('H','p','filename','ctrl'),
                    help="Create and run a regular wave.\nArguments:\n"
                    "H = wave height,\n"
                    "p = wave period,\n"
                    "filename = name of output file,\n"
                    "ctrl = control command.\n"
                    "Example:\n--regular_wave 1 6 output.csv linear "
                    )
parser.add_argument("--bretschneider_wave", nargs=4,metavar=('Hs','Te','filename','ctrl'),
                    help="Create and run a bretschneider sea state.\nArguments:\n"
                    "Hs = significant wave height\n"
                    "Te = energy period,\n"
                    "filename = name of output file,\n"
                    "ctrl = control command.\n"
                    "Example:\n--bretschneider_wave 1 6 output.csv linear "
                    )
parser.add_argument("--decay", nargs=4,metavar=('[x]','T','filename','ctrl'),
                    help="Run a decay test.\nArguments:\n"
                    "[x] = inital decay vector [z0,dz/dt0,x0,dx/dt0]\n"
                    "T = time of test,\n"
                    "filename = name of output file,\n"
                    "ctrl = control command.\n"
                    "Example:\n--decay '[1 0 1 0]' 20 'output.csv' 'linear' "
                    )
parser.add_argument("--GUI", action='store_true',
                    help="Start web server for GUI."
                    )
#args=parser.parse_args(["--regular_wave","3","3","lin","lin"])
args=parser.parse_args()
if args.regular_wave != None:
    a=args.regular_wave;
    reg_wave(float(a[0]),float(a[1]),a[2],a[3])
elif args.bretschneider_wave != None:
    a=args.bretschneider_wave;
    bretschneider_wave(float(a[0]),float(a[1]),a[2],a[3])
elif args.decay != None:
    a=args.decay;
    x=list(map(float,a[0][1:-1].replace(","," ").replace("  "," ").split()));
    decay_test(x,a[2],float(a[1]),a[3])
elif args.GUI:
    GUIServer.run()
#elif:
#print("Welcome to the COERbuoy1 model!\n")
#print("Usage:")
#print("regular wave: --regular_wave \n")
#print("2) Bretschneider sea state\n")
#print("3) Bretschneider sea state\n")
