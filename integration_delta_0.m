clc
clear
parameters_delta_0

dx=(L)/(Nx);
dy=H/(Ny);

% initial condition in terms of c(i,j) and ck(k)
for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j;
        c(i,j)=0;
        
        %initial condition for an impulsive injection
        c(3,j)=1/(Ny*dx*dy);
       

       ck(k)=c(i,j);
    end
end
ck=ck';

% setting up the integrator
dt=0.5;
tfinal=7;
tspan=[0:dt:tfinal];
Nt=floor(tfinal/dt);
ode_solver.options = odeset('reltol', 1e-3, 'abstol', 1e-6);
[t,ckn] = ode23(@vf_delta_0,tspan,ck,ode_solver.options);

%domain for plotting the solution
x=dx:dx:L;
y=dy:dy:H;
[X,Y]=meshgrid(x,y);

%collecting the concentration data in a 3D array
for l=1:Nt
    for i=1:Nx
        for j=1:Ny
            k=(i-1)*Ny+j;
            cn(i,j,l)=ckn(l,k);
        end
    end
end

% Plot Initial Condition Surface in a Separate Figure
figure; % Create new figure for initial condition
hold on;
surface(X', Y', c(:,:)); % Plot initial condition
shading interp; % Improve surface appearance
colorbar; % Add color scale
title('Initial Condition Surface');
xlabel('X');
ylabel('Y');
zlabel('c');
view(3); % Set 3D view
grid on;
hold off;

% Draw and update the simulation
figure; % Create new figure for time evolution
hold on;
for l = 2:2:Nt
    surface(X', Y', cn(:,:,l)); % Plot updated surface
    shading interp; % Smooth shading
    colorbar;
    title(['Time Evolution at Step ', num2str(l)]);
    xlabel('X');
    ylabel('Y');
    zlabel('c');
    view(3);
    grid on;
    drawnow; % Update figure dynamically
end
hold off;
