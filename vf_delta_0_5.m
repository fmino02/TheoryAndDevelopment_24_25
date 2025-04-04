function [B]=vf2d(t,ck)
parameters_delta_0.5

for i = 1:1:Nx
    for j = 1:1:Ny(i)
        k=sum(Ny(1:i-1))+j;
        a(k)=k;
    end
end

Lx=L/(Nx);
Ly=Lx*2*abs(delta)/L;
for j=1:Ny
    y(j)=j*Ly;
    %f(j)=0;
    Q(i,j) = 6 * H0 * (H2 * x - (L_y^2 * y(j) / 2) + yj * L_y - (L_y^3 * (H0 + delta * (2 * xi / L)) / 3) ...
    - (yj^2 * L_y * (H0 + delta * (2 * xi / L))) + (yj * L_y^2 * (H0 + delta * (2 * xi / L))));

end


%Volumes in contact with the inlet
for j=2:Ny-1
      vf(1,j)=-H/Lx/Ly*(c(1,j)*f(j)-c0*f(j)+H/Pe*(-(c(1,j+1)-c(1,j))*Lx/Ly+ ...
        (c(1,j)-c(1,j-1))*Lx/Ly-(c(2,j)-c(1,j))*Ly/Lx));
end

vf(1,1)=-H/Lx/Ly*(c(1,1)*f(1)-c0*f(1)+H/Pe*(-(c(1,2)-c(1,1))*Lx/Ly ...
    -(c(2,j)-c(1,j))*Ly/Lx));

vf(1,Ny)=-H/Lx/Ly*(c(1,Ny)*f(Ny)-c0*f(Ny)+H/Pe*((c(1,Ny)-c(1,Ny-1))*Lx/Ly- ...
    (c(2,Ny)-c(1,Ny))*Ly/Lx));


% internal internal volumes
for i=2:Nx-1
    for j=2:Ny-1
     vf(i,j)=-H/(Lx*Ly)*(c(i,j)*f(j)-c(i-1,j)*f(j)-H/Pe*((c(i,j+1)-2*c(i,j)+ ...
         c(i,j-1))*Lx/Ly+(c(i+1,j)-2*c(i,j)+c(i-1,j))*Ly/Lx));
    end
end

% Volumes in contact with the top boundary
j=Ny;
for i=2:Nx-1
    vf(i,j)=-H/(Lx*Ly)*(c(i,j)*f(j)-c(i-1,j)*f(j)+H/Pe*((c(i,j)-c(i,j-1))*Lx/Ly- ...
        (c(i+1,j)-c(i,j))*Ly/Lx+(c(i,j)-c(i-1,j))*Ly/Lx));
end
% Volumes in contact with the bottom boundary
j=1;
for i=2:Nx-1
    vf(i,j)=-H/(Lx*Ly)*(c(i,j)*f(j)-c(i-1,j)*f(j)+H/Pe*(-(c(i,j+1)-c(i,j))*Lx/Ly- ...
        (c(i+1,j)-c(i,j))*Ly/Lx+(c(i,j)-c(i-1,j))*Ly/Lx));
end
% volumes in contact with the outlet
i=Nx;
for j=2:Ny-1
    vf(i,j)=-H/Lx/Ly*(c(i,j)*f(j)-c(i-1,j)*f(j)+H/Pe*(-(c(i,j+1)-c(i,j))*Lx/Ly+ ...
        (c(i,j)-c(i,j-1))*Lx/Ly+(c(i,j)-c(i-1,j))*Ly/Lx));
end
% top right boundary volume
i=Nx;
j=Ny;
vf(i,j)=-H/Lx/Ly*(c(i,j)*f(j)-c(i-1,j)*f(j)+H/Pe*((c(i,j)-c(i,j-1))*Lx/Ly+(c(i,j)-c(i-1,j))*Ly/Lx));
% top left boundary volume
i=Nx;
j=1;
vf(i,j)=-H/Lx/Ly*(c(i,j)*f(j)-c(i-1,j)*f(j)+H/Pe*(-(c(i,j+1)-c(i,j))*Lx/Ly+(c(i,j)-c(i-1,j))*Ly/Lx));

for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j;
        B(k)=vf(i,j);
    end
end     
B=B';
     

     
     
     




       

    