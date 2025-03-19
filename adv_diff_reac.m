%A 1 d DVECTION_DIFFUSION_REACTION TRANSPORT EQUATION
clear all
close all

parametri_adr %read parameter values from a file
% 
% dx=(L)/(Nx+1); %length of the intervals in which we have divided the domain
% 
% % This is to integrate directly from 0 to 1
% C=zeros(1,Nx); % Initial Condition
% C(3)=1/dx; %impulsive concentration normalized so that the integral is 1
% tfinal=1;
% dt=0.08;
% Nt=floor(tfinal/dt); % number of time intervals
% tspan=[0:dt:tfinal]; % time instants at which we want the solution
% X=dx:dx:L-dx;
% %% First approach to the problem (saving information about a lot of time intervals
% % we use the integrator with its default values
%  [t,Cn] = ode23(@bilancio_adr,tspan,C);
% % [results of the integration as a matrix] = (vector field, time instants 
% % at which we want the solution, initial condition)
% % Cn is onethousand points at which we have discretized in space, 416 is
% % the amount of time instants considered
%  % plot(Cn(end,:),'-');
%  for i=2:Nt
%     plot(X',Cn(i,:),'-');
%     hold on;
%  end

%% Second approach to the problem
% define a custom ode solver specifying the parameters we need
 ode_solver.method = @ode23;
 ode_solver.options = odeset('reltol', 1e-3, 'abstol', 1e-6);
 ode_solve = @(tspan,C) ode_solver.method(@bilancio_adr, tspan, C, ode_solver.options);
% 
% 
% 
% Tmax=1;
% X=dx:dx:L-dx;
% dt=0.025;
% Nit=floor(Tmax/dt);
% told=0;
% % saving a lot of space because I am solving different problems with
% different initial conditions - storing at any time instant only the
% information we need. Saving a lot of space in terms of RAM, space needed
% to define the variables
% for jt=1:Nit  
%     [tnew,Cn]=ode_solve([told told+dt],C);
%     C=Cn(end,:);
%     told=tnew(end);
%     drawnow;
%     plot(X',C,'-');
%     hold on;
%    end