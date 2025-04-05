clc; clear; close all;

% Parameters
delta = 0.5;
H0    = 1;
L     = 10;
Nx    = 500;

% Generate Ny
Ny = Nx * H0 / delta + (1:Nx);              % Symmetric vertical resolution

% Grid spacing
dx = (L/2) / Nx;
dy = delta * dx / (L / 2);

% x-coordinates (cell boundaries)
x_flow=zeros(max(Ny),Nx);
    for j=1:Nx
    for i=1:Ny(j)
    x_flow(i,j)=dx*j;
    end
    end

% Generate y_grid
 % y coordinates stored in an array
 for i=1:Nx
     y_flow(1,i)=dy;
     for j=2:Ny(i)
         y_flow(j,i)=y_flow(j-1,i)+dy;
     end
 end

% Get H (channel height) per column
H = max(y_flow, [], 1);
%

% Interior left right flowrate (Q)
for i =1:1:Nx
    for j=1:1:Ny(i)
        Q(j,i) = 6 * H0 / (H(i))^2 * ( - dy^2 / 2 + y_flow(j,i) * dy ...
            - dy^3 / (3 * H(i)) - (y_flow(j,i)^2 * dy) / H(i) ...
            + y_flow(j,i) * dy^2 / H(i));
    end
end
for i = 1:Nx-1
    for j = 1:max(Ny)
        Qr(j,i) = Q(j,Nx-i);
    end
end
Q=[Q,Qr];
for i = 1:Nx-1
    for j = 1:max(Ny)
        x_flowr(j,i) = L - x_flow(j,Nx-i);
        if x_flowr(j,i) == L
            x_flowr(j,i) = 0;
        end
    end
end

x_flow=[x_flow,x_flowr];
for i = 1:Nx-1
    for j = 1:max(Ny)
        y_flowr(j,i) = y_flow(j,Nx-i);
    end
end

y_flow=[y_flow,y_flowr];

% Plot example horizontal flowrate profile
figure;
plot(y_flow(:,2),Q(:,2), 'LineWidth', 1.5)
xlabel('y')
ylabel('Horizontal flowrate')
title('Left-right flowrate at column i = ')
grid on

% Top-bottom flowrate evaluation
F_tb = zeros(max(Ny), Nx);  % Preallocate
%
j_eval = 1200;               % Chosen vertical index
for j = 1:Nx
    for i = 1:Ny(j)-1
        F_tb(i,j) = 6 * (y_flow(i,j))^2 * H0 * (-1 / (2 * (H(j))^2) ...
            + y_flow(i,j) / (3 * (H(j))^3) + 1 / (2 * (H(j)-2*delta*dx/L)^2) ...
            - y_flow(i,j) / (3 * (H(j)-2*delta*dx/L)^3));
    end
end
for i = 1:Nx-1
    for j = 1:max(Ny)
        F_tbr(j,i) = -F_tb(j,Nx-i);
    end
end

F_tb=[F_tb,F_tbr];

% Plot vertical flowrate at a specific height
figure;
plot(x_flow(j_eval,:), F_tb(j_eval,:), 'LineWidth', 1.5)
xlabel('x')
ylabel(['Vertical flowrate at y index j = ', num2str(j_eval)])
grid on