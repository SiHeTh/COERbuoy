## COERbuoy controller library

- **controller.m**, **controller.py**  
  Simple example for a COERbuoy controller written in python and octave  
  >python3 -m COERbuoy \--regular_wave 1 3 "file.csv" "controller.m"
  

- **MATLAB_controller.m**  
  Controller written in MATLAB. Starts the COERbuoy platform from within MATLAB to test the controller.  
  

- **controller_damping.py**  
  Applies the optimal damping for the selected regular wave 
  Example: optimal damping for a wave period of 6s:   
  >python3 -m COERbuoy \--regular_wave 1 6 "file.csv" "controller_damping.py 6s"
  
- **controller_reactive.py**
  Applies reactive control for the selected regular wave 
  Example: reactive control for a wave period of 6s:  
  >python3 -m COERbuoy  \--regular_wave 1 6 "file.csv" "controller_reactive.py 6s"
  

- **MomentBasedController.py**  
  Example of a linear moment based controller  
  Example: moment absed control with a stroke limit of +/- 1 m:  
  >python3 -m COERbuoy \--regular_wave 1 8 "file.csv" "MomentBasedController.py 1"
