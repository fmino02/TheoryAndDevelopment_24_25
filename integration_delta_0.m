clc
clear
parameters_delta_0

Nx
Nz
dx=(L)/(Nx);
dz=H/(Nz);



% initial condition in terms of c(i,j) and ck(k)
for i=1:Nx
    for j=1:Nz
        k=(i-1)*Nz+j;
        %c(i,j)=0;
        
        %initial condition for an impulsive injection
        c(3,j)=1/dx;
        ck(k)=c(i,j);
    end
end
ck=ck'

dt=1;
tfinal=10;
tspan=[0:dt:tfinal];
Nt=floor(tfinal/dt);
ode_solver.options = odeset('reltol', 1e-3, 'abstol', 1e-6);
[t,ckn] = ode23(@vf_delta_0,tspan,ck,ode_solver.options);
% ode_solver.method = @ode23;
% ode_solver.options = odeset('reltol', 1e-3, 'abstol', 1e-6);
% ode_solve = @(tspan,ckn) ode_solver.method(@vf_delta_0, tspan, ck, ode_solver.options);


x=dx:dx:L;
z=dz:dz:H;
[X,Y]=meshgrid(x,z);


for l=1:Nt
    for i=1:Nx
        for j=1:Nz
            k=(i-1)*Nz+j;
            cn(i,j,l)=ckn(l,k);
        end
    end
end

drawnow
view(3)
%plot initial condition
% surface(X',Y',c(:,:))
for l=2:2:Nt
  surface(X',Y',cn(:,:,l));
end
 hold off
