function [B] = vf_delta_negative(t, ck)
parameters_delta_negative;

% Fill concentration matrix
for i = 1:2*Nx
    for j = 1:Ny(i)
        k = sum(Ny(1:i-1)) + j;
        c(i,j) = ck(k);
    end
end

% y-coordinates
for i = 1:Nx
    y_flow(1,i) = dy;
    for j = 2:Ny(i)
        y_flow(j,i) = y_flow(j-1,i) + dy;
    end
end
for j=1:Ny(Nx)-1
y_flow(j,Nx+1)=dy*j;
end
H = max(y_flow, [], 1);  % Channel heights
Q=zeros(max(Ny),Nx);
% Compute Q
for i = 1:Nx
    for j = 1:Ny(i)-1
        Q(j,i) = 6 * H0 / H(i+1)^2 * ( -dy^2/2 + y_flow(j,i)*dy ...
            - dy^3/(3*H(i+1)) - y_flow(j,i)^2 * dy / H(i+1) ...
            + y_flow(j,i)*dy^2 / H(i+1));
    end
end

% for j = 1:Ny(Nx)-1
%         Q(j,Nx) = 6 * H0 / H(i+1)^2 * ( -dy^2/2 + y_flow(j,i)*dy ...
%             - dy^3/(3*H(i+1)) - y_flow(j,i)^2 * dy / H(i+1) ...
%             + y_flow(j,i)*dy^2 / H(i+1));
%     end

for i = 1:Nx-1
    for j = 1:max(Ny)-1
        Qr(j,i) = Q(j,Nx-i);
    end
    Qr(max(Ny),i)=0;
end

  for j = 1:max(Ny)
        Qr(j,Nx) = 6 / H0 * ( -dy^2/2 + y_flow(j,1)*dy ...
            - dy^3/(3*H0) - y_flow(j,1)^2 * dy / H0 ...
            + y_flow(j,1)*dy^2 / H0);
  end

Q = [Q, Qr];

Q = Q';

% Compute top-bottom flowrate
F_tb = zeros(max(Ny), Nx);
for j = 1:Nx
    for i = 1:Ny(j)-1
        F_tb(i,j) = 6 * y_flow(i,j)^2 * H0 * (-1/(2*H(j)^2) ...
            + y_flow(i,j)/(3*H(j)^3) + 1/(2*(H(j)-2*delta*dx/L)^2) ...
            - y_flow(i,j)/(3*(H(j)-2*delta*dx/L)^3));
    end
end
F_tbr =-flip(F_tb,2);
F_tb = [F_tb, F_tbr];

F_tb = F_tb';

vf = zeros(2*Nx, max(Ny));  % Preallocates full domain space
% Interior rectangular volumes
% First half of the domain
for i = 2:Nx-1
    for j = 2:Ny(i)
        vf(i,j) = -H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
             + c(i,j)*F_tb(i,j) - c(i,j+1)*F_tb(i,j+1) ...
             - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
             + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));
    end
end

% Second half of the domain
for i = Nx+1:2*Nx-1
    for j = 2:Ny(i)-1
        vf(i,j) = -H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
            + c(i,j)*F_tb(i,j) - c(i,j-1)*F_tb(i,j-1) ...
            - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
            + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));
    end
end

% Rectangular volumes at the inlet
i=1;
for j = 2:Ny(1)-1
        vf(i,j) = -H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
             + c(i,j)*F_tb(i,j) - c(i,j+1)*F_tb(i,j+1) ...
             - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
             + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));
end

% Rectangular volume in the bottom left corner
i=1;
j=1;
        vf(i,j) = -H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
             + c(i,j)*F_tb(i,j) - c(i,j+1)*F_tb(i,j+1) ...
             - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
             + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));

% Top triangles
%top left triangle
i = 1;
j = Ny(i);
        vf(i,j) = -2*H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
             + c(i,j)*F_tb(i,j) - c(i,j+1)*F_tb(i,j+1) ...
             - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
             + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));
%first half of triangles
for i = 2:Nx
    j = Ny(i);
        vf(i,j) = -2*H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
             + c(i,j)*F_tb(i,j) - c(i,j+1)*F_tb(i,j+1) ...
             - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
             + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));end
%second half of triangles
for i = Nx+1:2*Nx-1
    j = Ny(i);
        vf(i,j) = -2*H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
            + c(i,j)*F_tb(i,j) - c(i,j-1)*F_tb(i,j-1) ...
            - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
            + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));
end
%top right triangle
i = 2*Nx;
    j = Ny(i);
        vf(i,j) = -H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
            + c(i,j)*F_tb(i,j) - c(i,j-1)*F_tb(i,j-1) ...
            - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
            + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));

% Bottom rectangles
j = 1;
for i = 2:Nx
        vf(i,j) = -H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
             + c(i,j)*F_tb(i,j) - c(i,j+1)*F_tb(i,j+1) ...
             - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
             + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));
end

for i = Nx+1:2*Nx-1
        vf(i,j) = -H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
            + c(i,j)*F_tb(i,j) - c(i,j-1)*F_tb(i,j-1) ...
            - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
            + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));
end

% Outlet
i = 2*Nx;
for j = 2:Ny(i)-1
        vf(i,j) = -H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
            + c(i,j)*F_tb(i,j) - c(i,j-1)*F_tb(i,j-1) ...
            - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
            + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));
end
% rectangular volume in the bottom right corner
i=2*Nx;
j=1;
        vf(i,j) = -H0/(dx*dy)*(c(i,j)*Q(i,j)-c(i-1,j)*Q(i-1,j) ...
            + c(i,j)*F_tb(i,j) - c(i,j-1)*F_tb(i,j-1) ...
            - H0/Pe*((c(i,j+1)-2*c(i,j)+c(i,j-1))*dx/dy ...
            + (c(i+1,j)-2*c(i,j)+c(i-1,j))*dy/dx));
% Final output vector
B = [];
for i = 1:2*Nx
    for j = 1:Ny(i)
        k = sum(Ny(1:i-1)) + j;
        B(k,1) = vf(i,j);
    end
end
