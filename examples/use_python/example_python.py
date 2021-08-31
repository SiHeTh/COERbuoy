import numpy as np;
import COERbuoy;

##Examples with wave generator (beginner)
#regular wave with 1.5 m wave height and 3.5 s wave period; using the build-in linear damping
COERbuoy.reg_wave(1.5,3.5,"results/test.csv","linear");

#bretschneider wave with 1.5 m significant wave height and 12 s energy period; using the build-in linear damping
COERbuoy.bretschneider_wave(1.5,12,"results/bretschneider_wave.csv","linear");

#Regular wave with user defined controller
COERbuoy.reg_wave(3,12,"results/wave_with_controller.csv","python3 controller.py")#Linux/MacOS
COERbuoy.reg_wave(3,12,"results/wave_with_controller.csv","py controller.py")#Windows

#Regular wave, just open TCP (host) connection; the controller has to be run separately
#Example: Run the following line, then run "py controller.py"/"python3 controller.py" in a seperate process.
COERbuoy.reg_wave(3,12,"results/open_TCP.csv","TCP")

#Run a 30 s decay test with starting position 1 m above equilibrium position
COERbuoy.decay_test(1,"results/decay.csv",30,"linear")

#Run a period sweep
power=np.zeros(4);
for i,period in enumerate([4,6,8,12]):
    height=1;
    power[i]=COERbuoy.reg_wave(height,period,"results/regular_wave_H_1m_P_"+str(period)+"s.csv","linear");
print("Absorbed power for each sea state [kW]: "+str(power/1000));
print("Total absorbed power [kW]: "+str(np.sum(power/1000)));

##Calling start_simu directly (advanced)

#start_simu parameters (first is default):
#control = linear | TCP | <file command>       //controller to be used
#host    = True | False                        //TCP host or client
#name    = <filename>                          //name of the output file
#t0      = <seconds(float)>                    //transient time to wait until power measurement starting
#init    = [stroke,angular,v_stroke,v_angular] //initial conditions (float/SI units)
#EITHER
#file    = <csv-file>                          //file with wave data: time, wave-elevation
#OR 
#time    = <1xn list>                          //time series
#wave    = <1xn list>                          //corresponding wave elevation

#A user generated regular wave:
t = np.linspace(0,60,600);#time data
w = np.sin(0.1*t);#corresponding surface elevation
COERbuoy.start_simu(time=t, wave=w, name="results/user_generated.csv", t0=10, control="linear");

#Load a wave file
COERbuoy.start_simu(file="example_wave.csv", name="results/example_wave_output.csv", t0=10, control="linear");

