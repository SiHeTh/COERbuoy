clear -all
pkg load instrument-control
clc
 

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

gamma=100000; %[Ns/m] constant damping value
ctrl.F_pto = -gamma*status.stroke_speed(end); #set PTO force
ctrl.F_brake = 0; #set brake force

endfunction



%%-- do not change --%%

pause(4);
global BUFSIZE=sizeof(struct('time',zeros(100,1),'wave',zeros(100,1),'wave_forecast',zeros(100,1),'stroke_pos',zeros(100,1),'stroke_speed',zeros(100,1),'angular_pos',zeros(100,1),'angular_speed',zeros(100,1),'force',zeros(100,1),'test',zeros(100,1)));
global endian='B';
global s = tcpip('localhost', 5050, 'timeout', 100);
[type, mSize, endian] = computer
function [status]=receivemsg()
    global s;
    global endian;
    global BUFSIZE;
    array=tcp_read(s,BUFSIZE,10000);
    array=typecast(array,'double');
    if endian == 'L';
      array=swapbytes(array);
    endif
    array=array;
      
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
endfunction

%send control message
function sendmsg(time,force,brake,test)
    global endian;
    global s;
    array=double([time(1:9),force(1:9),brake(1:9),test(1:9)]);
    if endian == 'L'
      array=swapbytes(array);
    endif
    array=typecast(array,'uint8');
    tcp_write(s,uint8(array));
endfunction

%main
while 1
  #Constant damping control
  status=receivemsg();
  [ctrl]=control(status);
  sendmsg(linspace(status.time(end),status.time(end)+1,9),ones(1,9)*ctrl.F_pto,ones(1,9)*ctrl.F_brake,zeros(1,9));
endwhile


