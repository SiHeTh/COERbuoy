Example of a custom WEC for the COERbuoy platform

1. How it works
2. The settings file
3. The WEC data folder

1. How it works
When started, the COERbuoy platform searches the working directory, usually the folder the program was started from, for a file called "coerbuoy_settings.txt". If this file exists, it will use the WEC model and settings provided in this file instead of using the default configuration.
To run the example, make sure to start the COERbuoy platform from within this folder. This can be done by open the command line / terminal and navigate to this folder. Then type COERbuoy.GUIServer to start the COERbuoy platform. Open a web browser and open "localhost:8080". Navigate to "Parameters". The parameters text box to the left should show the parameters of the custom WEC. The example will furthermore show a text highlighting that this is a custom WEC.

2. The settings file:
The setting files is a json formatted text file.
"host" specifies the address for the GUI server (default: "localhost")
"port" is the port number for the GUI server (default: 8080)
"hydro" specifies the hydrodynamic model (default: "Floater_BEM")
"WECfolder" folder that specifies the WEC (see chapter 3)
"resolution" time resolution (in seconds) for the results file (default: 0.01)
"status_message" array specifying which elements of the status message are send (default: [1,1,1,1,1,1,1,1,1,1])

3. The WEC data folder
The structure of the WEC data folder is:

WEC_data
|
|-floater.txt #WEC parameters
|-dynamics.py #WEC logic
|
|-BEM #linear hydrodynamic data
    |h_-1_p_0 #... 1 m submerged, no pitch
        |....
    |h_0_p_0 #... heave mean position, no pitch
        |....
    |h_1_p_0 #... -1 m submerged, no pitch
        |
        |-results
            |
            |DiffractionForce.tec #nemoh standard
            |ExcitationForce.tec #nemoh standard
            |FKforce.tec #nemoh standard
            |RadiationCoefficents.tec #nemoh standard
            |IRF.tec #nemoh standard

3.1 dynamics.py
This file describes the logic of the WEC. See the file for further information.

3.2 floater.txt
Specifies the parameter used by dynamics.py. For the calculation of the hydrostatic force, it must specify the geometry in terms of truncated cone section, where each section is set by four coordinates ("coord":[x1,r1,x2,r2]):
x1: the upper end (relative to the mean water level) of the cone section
r1: the radius at x1
r2: the lower end (relative to the mean water level) of the cone section
r2: the radius at x2

The first section must start with r2=0; the last section must end with r1=0; For two consecutive sections i and j the radii r1_i=r2_j and x coordinates x1_i=x2_j must match.

Example:
"geo":[{"type":"cone","coord":[-4.5,0,-2,3]},
{"type":"cone","coord":[-2,3,-1,3.7]},
{"type":"cone","coord":[-1,3.7,0,4]},
{"type":"cone","coord":[0,4,1,3.7]},
{"type":"cone","coord":[1,3.7,2,3]},
{"type":"cone","coord":[2,3,4.5,0]}],


3.3 BEM
The files here follow the NemoH standard. While the number of frequencies is variable, all files must be calculated for three degrees of freedom (heave, surge, pitch). Four lines of header are assumed. By now, Floater_BEM can only handle heave and surge, thus, all files must be calculated for zero pitch. The files must be evenly spaced in terms of submergence. For further details, see the example data.
