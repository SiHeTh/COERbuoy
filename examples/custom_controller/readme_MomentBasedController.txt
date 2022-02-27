Basic moment domain controller

This is a self-parametrising, one degree-of-freedom, linear receeding-horizon moment based controller for the COERbuoy platform with stroke limit, see [1].

To use it in COERbuoy there are several option:

- Look for the folder "COERbuoy_data/controller" in your home directory and
  copy the file their, then start the GUI and select
  "MomentBasedController.py".
- open a terminal, navigate to the directory where the control file is 
  located and run COERbuoy from command line, using
  "MomentBasedController.py" as controller argument.
- start the COERbuoy simulation with the control command "TCP", which opens
  the ControlInterface, the run "MomentBasedController.py".


[1] N. Faedo et al.: Energy-maximising control of wave energy converters
using a moment-domain representation, Control Engineering Practice 81, 
2018
