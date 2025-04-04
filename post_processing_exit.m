%% Post processing and comparison of MATLAB and COMSOL simulations
clear
clc

% === Load MATLAB simulation data ===
filename_matlab = "delta_0_Pe_100_matlab.xlsx";
data_matlab = readtable(filename_matlab);
c_matlab = data_matlab.right_concentration;
t_matlab = data_matlab.Time;

% === Load COMSOL simulation data ===
filename_comsol = "delta_0_Pe_100_comsol.xlsx";
data_comsol = readtable(filename_comsol);
c_comsol = data_comsol.delta_0_model_mph;
t_comsol = data_comsol.Model;

% === Normalize both datasets to M_in = 1 ===
M_in = 1;

% Process MATLAB data
l_matlab = 50;
M_out_matlab = trapz(t_matlab(1:l_matlab), c_matlab(1:l_matlab));
while M_out_matlab < M_in && l_matlab < length(t_matlab)
    l_matlab = l_matlab + 1;
    M_out_matlab = trapz(t_matlab(1:l_matlab), c_matlab(1:l_matlab));
end
t_sel_matlab = t_matlab(1:l_matlab);
c_sel_matlab = c_matlab(1:l_matlab);
m1_matlab = trapz(t_sel_matlab, c_sel_matlab .* t_sel_matlab);
m2_matlab = trapz(t_sel_matlab, c_sel_matlab .* t_sel_matlab.^2);
m3_matlab = trapz(t_sel_matlab, c_sel_matlab .* t_sel_matlab.^3);
var_matlab = m2_matlab - m1_matlab^2;
skw_matlab = (m3_matlab - 3 * m1_matlab * var_matlab - m1_matlab^3) / (var_matlab^(3/2));

% Process COMSOL data
l_comsol = 50;
M_out_comsol = trapz(t_comsol(1:l_comsol), c_comsol(1:l_comsol));
while M_out_comsol < M_in && l_comsol < length(t_comsol)
    l_comsol = l_comsol + 1;
    M_out_comsol = trapz(t_comsol(1:l_comsol), c_comsol(1:l_comsol));
end
t_sel_comsol = t_comsol(1:l_comsol);
c_sel_comsol = c_comsol(1:l_comsol);
m1_comsol = trapz(t_sel_comsol, c_sel_comsol .* t_sel_comsol);
m2_comsol = trapz(t_sel_comsol, c_sel_comsol .* t_sel_comsol.^2);
m3_comsol = trapz(t_sel_comsol, c_sel_comsol .* t_sel_comsol.^3);
var_comsol = m2_comsol - m1_comsol^2;
skw_comsol = (m3_comsol - 3 * m1_comsol * var_comsol - m1_comsol^3) / (var_comsol^(3/2));

% === Plot both chromatograms ===
figure;
plot(t_matlab, c_matlab, 'b', 'LineWidth', 2);
hold on;
plot(t_comsol, c_comsol, 'r', 'LineWidth', 2);
xlabel('Tempo');
ylabel('Concentrazione media outlet');
title('Confronto Chromatogrammi: MATLAB vs COMSOL');
legend('MATLAB', 'COMSOL', 'Location', 'Best');
grid on;

% === Save the moment results to a table ===
% comparisonTable = table( ...
%     ["MATLAB"; "COMSOL"], ...
%     [m1_matlab; m1_comsol], ...
%     [var_matlab; var_comsol], ...
%     [
