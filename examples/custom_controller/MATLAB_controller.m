clear -all
clc
%% Autorun controller for MATLAB for the COERbuoy platform and the COERbuoy1 benchmark
% 2021 Centre for Ocean Energy Research, Maynooth University
%
% This controller is a template for using the COERbuoy ControlInterface
% with MATLAB. For demonstartion it implements a simple linear damping
% controller.
% Requirements: 1) python3 and the COERbuoy package have to be installed
%               2) MATLAB > 2016 needs to be installed
%               3) Set the correct command for python in (A)
%
% What it does: 1) COERbuoy is run in a terminal with the settings specified in (A)
%               2) Establishes a connection to COERbuoy in (B)
%               3) Then starts the main loop in (C)
%               4) And calls the control algorithm (D)
%
% How to implement controllers: Change the controller logic in (D).
% For questions and ideas, please write to simon.thomas.2021@mumail.ie

%% (A) Start COERbuoy simulation with desired parameters

%Usage of StartCOERbuoy
%StartCOERbuoy(cmd, spectra, height, period, filename)
    % cmd      - python interpreter command - python3 | py | python
    % spectra  - wave spectrum              - regular | irregular
    % height   - (significant) wave height  - 0.1 - 10.0
    % period   - wave (energy) period       - 0.1 - 20.0
    % filename - output filename            - string
    
% TODO: Check the cmd field corresponds with your system
% Deactivate this function when you call the controller from the COERbuoy platform
StartCOERbuoy("python3", "regular" , 1, 8, "output.csv")


%% (B) Initialize();
Initialize();


%% (C) Main loop
while 1
  %Constant damping control
  status=receivemsg();
  [ctrl]=control(status);
  sendmsg(linspace(status.time(end),status.time(end)+1,9),ones(1,9)*ctrl.F_pto,ones(1,9)*ctrl.F_brake,zeros(1,9));
end

%% (D) The control algorithm - Please replace with your own ideas
function [ctrl]=control(status)
    %status.time          - a sequence of 100 values, the last element contains the most recent value
    %status.wave          - a sequence of 100 values with the wave elevation corresponding to the time series*
    %status.wave_forecast - a sequence of 100 values with the wave forecast*
    %status.stroke_pos    - a sequence of 100 values with the stroke position; correspoding to the time series
    %status.stroke_speed  - a sequence of 100 values with the stroke speed; corresponding to the time series
    %status.angular_pos   - a sequence of 100 values with the angular position; corresponding to the time series
    %status.angular_speed - a sequence of 100 values with the angular speed; corresponding to the time series
    %status.force         - a sequence of 100 values with the force measured by the force sensor; corresponding to the time series
    %* not available in the COERbuoy1 benchmark

    %Example of a linear damping control
    gamma=300000; %[Ns/m] constant damping value
    ctrl.F_pto = -gamma*status.stroke_speed(end); %set PTO force
    ctrl.F_brake = 0; %set brake force

    end



%% Controlinterface connection handling - DO NOT CHANGE
% DO NOT CHANGE -- DO NOT CHANGE -- DO NOT CHANGE -- DO NOT CHANGE

%Start the COERbuoy platform
function StartCOERbuoy(py,wave,H,P,outfile)%python_cmd, regular | irregular, (significant) wave height, wave (energy) period, filename
    wavecmd="--bretschneider_wave ";
    if strcmp(wave,"regular")
        wavecmd="--regular_wave ";
    end
    wavecmd=wavecmd+string(H)+" "+string(P);
    
    cmd=py+" -m COERbuoy "+wavecmd+' "'+outfile+'" "TCP" &';
    %Start COERbuoy
    system(cmd)
end

%Initialize the ControlInterface TCP connection
function Initialize()
    global BUFSIZE;
    global endian;
    global s;
    pause(2)
    BUFSIZE=100*8*9;
    endian='B';
    s = tcpclient('localhost', 5050);
    [type, mSize, endian] = computer;
end

%send control message
function sendmsg(time,force,brake,test)
    global endian;
    global s;
    array=double([time(1:9),force(1:9),brake(1:9),test(1:9)]);
    if endian == 'L'
      array=swapbytes(array);
    end
    array=typecast(array,'uint8');
    write(s,uint8(array));
end
    
%receive control message
function [status]=receivemsg()
    global s;
    global endian;
    global BUFSIZE;
    array=read(s,BUFSIZE);
    array=typecast(array,'double');
    if endian == 'L'
      array=swapbytes(array);
    end
      
    i=1;
    status.time=array(1+(i-1)*100:i*100);
    i=i+1;
    status.wave=array(1+(i-1)*100:i*100);
    i=i+1;
    status.wave_forecast=array(1+(i-1)*100:i*100);
    i=i+1;
    status.stroke_pos=array(1+(i-1)*100:i*100);
    i=i+1;
    status.stroke_speed=array(1+(i-1)*100:i*100);
    i=i+1;
    status.angular_pos=array(1+(i-1)*100:i*100);
    i=i+1;
    status.angular_speed=array(1+(i-1)*100:i*100);
    i=i+1;
    status.force=array(1+(i-1)*100:i*100);
    i=i+1;
    status.test=array(1+(i-1)*100:i*100);
end


