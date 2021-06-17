from .simulation import start_simu, reg_wave, bretschneider_wave, decay_test
import argparse
from argparse import RawTextHelpFormatter as tf;

parser=argparse.ArgumentParser(description="The COERbuoy model",formatter_class=tf)
parser.add_argument("--regular_wave", nargs=4,type=str, metavar=('H','p','filename','ctrl'),
                    help="Create and run a regular wave.\nArguments:\n"
                    "H = wave height,\n"
                    "p = wave period,\n"
                    "filename = name of output file,\n"
                    "ctr = control command.\n"
                    "Example:\n--regular_wave 1 6 output.csv linear "
                    )
parser.add_argument("--bretschneider_wave", nargs=4,metavar=('Hs','Te','filename','ctrl'),
                    help="Create and run a bretschneider sea state.\nArguments:\n"
                    "Hs = significant wave height\n"
                    "Te = energy period,\n"
                    "filename = name of output file,\n"
                    "ctr = control command.\n"
                    "Example:\n--bretschneider_wave 1 6 output.csv linear "
                    )
                    
#args=parser.parse_args(["--regular_wave","3","3","lin","lin"])
args=parser.parse_args()
print(args)
if args.regular_wave != None:
    a=args.regular_wave;
    reg_wave(float(a[0]),float(a[1]),a[2],a[3])
elif args.bretschneider_wave != None:
    a=args.bretschneider_wave;
    bretschneider_wave(float(a[0]),float(a[1]),a[2],a[3])
#elif:
#print("Welcome to the COERbuoy1 model!\n")
#print("Usage:")
#print("regular wave: --regular_wave \n")
#print("2) Bretschneider sea state\n")
#print("3) Bretschneider sea state\n")
