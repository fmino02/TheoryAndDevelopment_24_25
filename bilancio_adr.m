function [VF] = bilancio_adr(t,C) %this function is reporting a vector field called VF
% vector fields need to be written as nonautonomous because the ODE
% integrator wants them so
parametri_adr %read parameter values from a file
dx=1/(Nx+1);
% internal internal nodes - we are considering the nonlinear problem
for i=2:Nx-1
    VF(i)=(C(i+1)-2*C(i)+C(i-1))/dx^2/Pe-(C(i)-C(i-1))/dx-Da*C(i)^n;
end
% vector field at the boundary nodes
VF(1)=(C(2)-2*C(1)+Cin)/dx^2/Pe-(C(1)-Cin)/dx-Da*C(1)^n;
VF(Nx)=(-C(Nx)+C(Nx-1))/dx^2/Pe-(C(Nx)-C(Nx-1))/dx-Da*C(Nx)^n;
VF=VF'; %trick to have a column vector, which the integrator wants
end

