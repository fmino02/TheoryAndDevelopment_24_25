%A 2 d DIFFUSION_REACTION TRANSPORT EQUATION
clear all
close all

Nx=50;
Ny=50;
Lx=1;
Ly=1;
Qtilde=2;
Tin=1;
Tout=0;

dx=Lx/Nx;
dy=Ly/Ny;

% INITIAL CONDITION FOR T(i,j)
for i=1:Nx
    for j=1:Ny
        MT(i,j)=Tout;
    end
end
% INITIAL CONDITION FOR T(k)
for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j
        VT(k)=MT(i,j);
    end
end



% This is to integrate directly from 0 to 1
tfinal=2;
dt=0.2;
Nt=floor(tfinal/dt);
tspan=[0:dt:tfinal];

[t,VTn] = ode23(@balance2d,tspan,VT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=dx/2:dx:Lx-dx/2
y=dy/2:dy:Ly-dy/2
[X,Y]=meshgrid(x,y)


for l=1:Nt
    for i=1:Nx
        for j=1:Ny
            k=(i-1)*Ny+j;
            MTp(i,j,l)=VTn(l,k);
        end
    end
end

drawnow
view(3)
 surface(X',Y',MTp(:,:,2))
%  hold off

