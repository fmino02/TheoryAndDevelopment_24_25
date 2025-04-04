%% Flowrate for delta > 0
H0 = 1;
L = 10 * H0;
Nx = 50;
Ny = 30;
Lx = L / Nx;
Ly = H0 / Ny;

% Preallocate arrays
y = zeros(1, Ny);
Q = zeros(1, Ny);

% Compute y values
for j = 1:Ny
    y(j) = j * Ly;
end

% Compute Q values and plot
figure;
hold on;
for j = 1:Ny
    Q(j) = (6 / H0) * (- (Ly^2) / 2 + y(j) * Ly - (Ly^3) / (3 * H0) ...
        - (y(j)^2 * Ly) / H0 + (y(j) * Ly^2) / H0);
    
    plot(y(j), Q(j), 'ro'); % Plot each point as a red circle
end
hold off;
xlabel('y');
ylabel('Q');
title('Flowrate Distribution');
grid on;
