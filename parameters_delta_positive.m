clc
H0=1; %scaling height - H0
L=10*H0; %length (axial) of the channel
c0=0; %inlet concentration
Pe=10; %need to study the interval (1-1000)
Nx=5; %number of divisions in the x direction
delta=0.5*H0;
dx=(L/2)/(Nx);
%per rispettare il rapporto tra i lati dei cateti il dy del triangolo è
%fissato
dy=dx*2*delta/L;
% Generate Ny
Ny = Nx * H0 / delta + (1:Nx);
Ny=[Ny,flip(Ny)];


