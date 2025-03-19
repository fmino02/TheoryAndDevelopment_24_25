%DIFFUSION IN a PORE
parametri

dx=1/(Nx-1);
dy=1/(Ny-1);

% INITIAL CONDITION 
C=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j;
        CV(k)=C(i,j);
    end
end

ode_solver.method = @ode23;
ode_solver.options = odeset('reltol', 1e-3, 'abstol', 1e-6);
ode_solve = @(tspan,CV) ode_solver.method(@bilancio_pore, tspan, CV, ode_solver.options);


Tmax=0.3;
dt=0.01;
Nit=fix(Tmax/dt);
told=0;
for jt=1:Nit
    [tnew,CVn]=ode_solve([told told+dt],CV);
    CV=CVn(end,:);
    told=tnew(end)
    told
    for i=1:Nx
        for j=1:Ny
            k=(i-1)*Ny+j;
            CC(i+1,j)=CV(k);
        end
    end
    for j=1:Ny
        CC(1,j)=1;
    end
    drawnow;
   
    %axis([0 250 0 50 0 1]);
    surfc(CC,'FaceColor','interp','EdgeColor','none');
    %contourf(CC);
    hold off;
end
