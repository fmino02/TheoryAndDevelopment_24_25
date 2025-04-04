clc
clear

% Parameters
delta = 1;
H0 = 1;
L = 10;
Nx0 = 25;

% Compute Nx
Nx = 2*Nx0 * delta / H0;
Nx_full = 2 * Nx;

% Define Ny symmetric profile
Ny = zeros(1, Nx_full);
for i = 1:Nx
    Ny(i) = 2*Nx0 + i;
    Ny(Nx_full - i + 1) = 2*Nx0 + i;
end
Nx = Nx_full;

% Grid spacing
dx = L / Nx;
dy = delta * dx / (L / 2);
%% Coordinates at which concentrations should be evaluated
% x coordinates (cell centers)
x_conc = dx/2 : dx : L - dx/2;

% y coordinates stored in an array
for i=1:Nx
    y_conc(1,i)=dy/2;
    for j=2:Ny(i)
        y_conc(j,i)=y_conc(j-1,i)+dy;
    end
end
%% Coordinates at which the flowrate should be evaluated
% x coordinates (boundaries of the cells)
x_flow = 0 : dx : L;

 % y coordinates stored in an array
 for i=1:Nx
     y_flow(1,i)=0;
     for j=2:Ny(i)+1
         y_flow(j,i)=y_flow(j-1,i)+dy;
     end
 end

 H = max(y_flow, [], 1);

for i =1:1:Nx
    for j=1:1:Ny(i)
        Q(i,j) = 6 * H0 / ((H(i))^2 * ( - dy^2 / 2 + y_flow(i,j) * dy ...
            - dy^3 / (3 * H(i)) - (y_flow(i,j)^2 * dy) / H(i) ...
            + (y_flow(i,j) * dy^2) / H(i)));
    end
end
A=x_flow(1:end-1);
B=Q(:,3);
plot(A,B)


