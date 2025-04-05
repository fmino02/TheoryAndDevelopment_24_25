clc
clear

% Parameters
delta = -0.5;
H0 = 1;
L = 10;
Nx0 = 100;

% Compute Nx
Nx = 2*Nx0 * abs(delta) / H0;
Nx_full = 2 * Nx;

% Define Ny symmetric profile
Ny = zeros(1, Nx_full);
for i = 1:Nx
    Ny(i) = 2*Nx0  - i;
    Ny(Nx_full - i) = 2*Nx0 - i;
end
Nx = Nx_full;

% Grid spacing
dx = L / Nx;
dy = abs(delta) * dx / (L / 2);
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
x_flow = dx : dx : L;

 % y coordinates stored in an array
 for i=1:Nx
     y_flow(1,i)=dy;
     for j=2:Ny(i)+1
         y_flow(j,i)=y_flow(j-1,i)+dy;
     end
 end

 H = max(y_flow, [], 1);

y__flow=y_flow';

for i =1:1:Nx
    for j=1:1:Ny(i)+1
        Q(i,j) = 6 * H0 / (H(i))^2 * ( - dy^2 / 2 + y__flow(i,j) * dy ...
            - dy^3 / (3 * H(i)) - (y__flow(i,j)^2 * dy) / H(i) ...
            + y__flow(i,j) * dy^2 / H(i));
    end
end


    for j=1:1:Ny(i)+1
        F(1,j) = 6 * (y__flow(1,j))^2 * H0 * (-1 / (2 * (H(1))^2) ...
            + y__flow(1,j) / (3 * (H(1))^3) + 1 / (2 * (H0)^2) ...
            - y__flow(i,j) / (3 * (H0)^3));
    end


for i =2:1:Nx
    for j=1:1:Ny(i)+1
        F(i,j) = 6 * (y__flow(i,j))^2 * H0 * (-1 / (2 * (H(i))^2) ...
            + y__flow(i,j) / (3 * (H(i))^3) + 1 / (2 * (H(i-1))^2) ...
            - y__flow(i,j) / (3 * (H(i-1))^3));
    end
end


figure;
A=x_flow(1:(end));
B=Q(:,3);
plot(A,B)
xlabel('x');
ylabel('Portata al bordo destro');
grid on;

figure;
C=y_flow(1:(end),2);
C=C';
D=F(3,:);
G=F(:,3);
plot(C,D)
xlabel('y');
ylabel('Portata al bordo in alto');
grid on;

figure
plot(A,G)
xlabel('x');
ylabel('Portata al bordo in alto');
grid on;


