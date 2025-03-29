function [B]=vf2d(t,ck)
parameters_delta_0

for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j;
        c(i,j)=ck(k);
    end
end

Lx=L/(Nx);
Ly=H/(Ny);
for j=1:Ny
    y(j)=j*Ly;
    f(j) = (6 / H) * (- (Ly^2) / 2 + y(j) * Ly - (Ly^3) / (3 * H) ...
        - (y(j)^2 * Ly) / H + (y(j) * Ly^2) / H);
end

% internal internal volumes
for i=2:Nx-1
    for j=2:Ny-1
     vf(i,j)=-H/(Lx*Ly)*(-c(i,j)*f(j)-c(i-1,j)*f(j)+H/Pe*((c(i,j+1)-2*c(i,j)+ ...
         c(i,j-1))*Lx/Ly+(c(i+1,j)-2*c(i,j)+c(i-1,j)*Ly/Lx)));
    end
end
% Volumes in contact with the inlet
i=1;
for j=1:Ny
    c(i,j)=c0;
end
% Volumes in contact with the top boundary
j=Ny;
for i=2:Nx-1
    vf(i,j)=-H/Lx/Ly*(c(i,j)*f(j)+c(i-1,j)*f(j)-H/Pe*((c(i,j)-c(i,j-1))*Lx/Ly- ...
        (c(i+1,j)-c(i,j))*Ly/Lx+(c(i,j)-c(i,j-1))*Lx/Ly+(c(i,j)-c(i-1,j))*Ly/Lx));
end
% Volumes in contact with the bottom boundary
j=1;
for i=2:Nx-1
    vf(i,j)=-H/(Lx*Ly)*(c(i,j)*f(j)-c(i-1,j)*f(j)-H/Pe*((c(i,j+1)-c(i,j))*Lx/Ly- ...
        (c(i+1,j)-c(i,j))*Ly/Lx-(c(i,j)-c(i-1,j))*Ly/Lx));
end
% volumes in contact with the outlet
i=Nx;
for j=2:Ny-1
    vf(i,j)=-H/Lx/Ly*(c(i,j)*f(j)-c(i-1,j)*f(j)-H/Pe*((c(i,j+1)-c(i,j))*Lx/Ly- ...
        (c(i,j)-c(i,j-1))*Lx/Ly-(c(i,j)-c(i-1,j))*Ly/Lx));
end
% boundary point E
i=Nx;
j=Ny;
vf(i,j)=-H/Lx/Ly*(c(i,j)*f(j)-c(i-1,j)*f(j)-H/Pe*(-(c(i,j)-c(i,j-1))*Lx/Ly-(c(i,j)-c(i-1,j))*Ly/Lx));
% boundary point F
i=Nx;
j=1;
vf(i,j)=-H/Lx/Ly*(c(i,j)*f(j)-c(i-1,j)*f(j)-H/Pe*((c(i,j+1)-c(i,j))*Lx/Ly-(c(i,j)-c(i-1,j))*Ly/Lx));

for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j;
        B(k)=vf(i,j);
    end
end     
B=B';
     

     
     
     




       

    