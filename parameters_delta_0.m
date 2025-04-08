H=1; %scaling height - H0
L=10*H; %length (axial) of the channel
c0=0; %inlet concentration
Pe=1; %need to study the interval (1-1000)
Nx=80; %number of divisions in the x direction
Ny=20; %number of divisions in the y direction
start=floor(Nx*0.1);
if start==1
    start=start+1;
end
