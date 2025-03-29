clc
clear
parameters_delta_0

Nx
Ny
dx=(L)/(Nx);
dy=H/(Ny);



% initial condition in terms of c(i,j) and ck(k)
for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j;
        c(i,j)=0;
        
        %initial condition for an impulsive injection
        c(3,j)=1/(H*dx);
       

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
y=dy:dy:H;
[X,Y]=meshgrid(x,y);


for l=1:Nt
    for i=1:Nx
        for j=1:Ny
            k=(i-1)*Ny+j;
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
