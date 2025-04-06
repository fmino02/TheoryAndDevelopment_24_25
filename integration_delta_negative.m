clc
clear
close all
parameters_delta_negative

ck = zeros(sum(Ny), 1);

% initial condition in terms of c(i,j) and ck(k)
for i = 1:1:2*Nx
    for j = 1:1:Ny(i)
        k=sum(Ny(1:i-1))+j;
        c(i,j)=0;
        
        %initial condition for an impulsive injection
        c(3,j)=1/(Ny(3)*dx*dy);
       

        ck(k)=c(i,j);
    end
end
ck=ck';

% setting up the integrator
dt=0.05;
tfinal=50;
tspan=[0:dt:tfinal];
Nt=floor(tfinal/dt);
ode_solver.options = odeset('reltol', 1e-3, 'abstol', 1e-6);
[t,ckn] = ode23(@vf_delta_negative,tspan,ck,ode_solver.options);

%domain for plotting the solution
x=dx/2:dx:L-dx/2;
y=dy/2:dy:max(Ny)*dy;
[X,Y]=meshgrid(x,y);

%collecting the concentration data in a 3D array
for l=1:Nt
    for i = 1:1:2*Nx
        for j = 1:1:Ny(i)
        k=sum(Ny(1:i-1))+j;
            cn(i,j,l)=ckn(l,k);
        end
    end
    H(l)=sum(sum(cn(:,:,l)).*dx.*dy);
end

for l=1:Nt
    for j=1:Ny(2*Nx)
        C(j,l)=cn(2*Nx,j,l);
    end
    D(l)=sum(C(:,l))/Ny(2*Nx);
    T(l)=l*dt;
end

% Save to table for post-processing
resultsTable = table(T', D', 'VariableNames', {'Time', 'right_concentration'});
writetable(resultsTable, 'delta_0_Pe_100_matlab.xlsx');

D1=sum(D.*(1:Nt)*dt)*dt;
D2=sum(D.*(1:Nt).^2*dt*dt)*dt;
D3=sum(D.*(1:Nt).^3*dt*dt*dt)*dt;
Var=D2-D1^2;
Kur=(D3-3*D1*Var-D1^3)/(Var^(3/2));

figure;
plot(T, D, 'b', 'LineWidth', 2);
xlabel('Tempo');
ylabel('Concentrazione media outlet');
title('Grafico di T vs. D');
grid on;

figure;
plot(T, H, 'r', 'LineWidth', 2);
xlabel('Tempo');
ylabel('Massa totale');
title('Grafico di T vs. H');
grid on;

% % Plot Initial Condition Surface in a Separate Figure
% figure; % Create new figure for initial condition
% hold on;
% surface(X', Y', c(:,:)); % Plot initial condition
% shading interp; % Improve surface appearance
% colorbar; % Add color scale
% title('Initial Condition Surface');
% xlabel('X');
% ylabel('Y');
% zlabel('c');
% view(3); % Set 3D view
% grid on;
% hold off;

% % Draw and update the simulation
% figure; % Create new figure for time evolution
% hold on;
% for l = 2:2:Nt
%     surface(X', Y', cn(:,:,l)); % Plot updated surface
%     %shading interp; % Smooth shading
%     %colorbar;
%     title(['Time Evolution at Step ', num2str(l)]);
%     xlabel('X');
%     ylabel('Y');
%     zlabel('c');
%     view(3);
%     grid on;
%     drawnow; % Update figure dynamically
% end
% hold off;
