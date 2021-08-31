COERbuoy custom controller (python/octave)

1. How it works
2. Debugging the controller

1. How it works

1.1 Controller written in python
Start python from the command line with the controller as argument:
  python3 -m COERbuoy --regular_wave 1 2 "test.csv" "python3 controller.py"
  or
  py -m COERbuoy --regular_wave 1 2 "test.csv" "py controller.py"

The COERbuoy has to be started from within the same folder as "controller.py". Alternatively, the full path can be used:
  python3 -m COERbuoy --regular_wave 1 2 "test.csv" "python3 /home/user/WEC/controller.py"
  or
  py -m COERbuoy --regular_wave 1 2 "test.csv" "py C:\Users\Me\WEC\controller.py"


1.2 Controller written in octave (and similar MATLAB)
Start python from the command line with the controller as argument:
  python3 -m COERbuoy --regular_wave 1 2 "test.csv" "octave controller.m"
  or
  py -m COERbuoy --regular_wave 1 2 "test.csv" "octave controller.m"

The COERboy has to be started from within the same folder as "controller.py". Alternatively, the full path can be used:
  python3 -m COERbuoy --regular_wave 1 2 "test.csv" "octave /home/user/WEC/controller.m"
  or
  py -m COERbuoy --regular_wave 1 2 "test.csv" "octave C:\Users\Me\WEC\controller.m"


2. Debugging the controller
For convenience, the COERbuoy platform will handle the TCP connection, stating the controller, connect and start the data exchange, automatically. However, when debugging, it may be useful to run the controller in a debugger. One way to achieve this, is to start COERbuoy in TCP mode, so that it only opens a TCP connection, but does not run the controller. The controller can then be run externally, for example inside an IDE.
1. Run the controller using "TCP" instead of the controller name:
   >py -m COERbuoy --regular_wave 1 2 "test.csv" "TCP"
   COERboy will open  TCP connection and wait for 10 seconds that a controller connects.
2. Start the controller from the command line / IDE.
3. Both should connect and COERbuoy should show "Start solver"
   If an error 'timed out' occur, the connection was closed before the controller was connected. This can either indicate a problem with the controller, or the controller was started too late.


