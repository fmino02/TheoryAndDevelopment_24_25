function [VF]=balance2d(t,VT)

Nx=50;
Ny=50;
Lx=1;
Ly=1;
Qtilde=2;
Tin=1;
Tout=0;

dx=Lx/Nx;
dy=Ly/Ny;

for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j;
        MT(i,j)=VT(k);
    end
end

% A REGION
for i=2:Nx-1
    for j=2:Ny-1
        MF(i,j)=(-(MT(i,j)-MT(i-1,j))/dx^2)-(-(MT(i+1,j)-MT(i,j))/dx^2)+...
            (-(MT(i,j)-MT(i,j-1))/dy^2)-(-(MT(i,j+1)-MT(i,j))/dy^2);
    end
end
% B REGION
i=1;
for j=2:Ny-1
    MF(i,j)=(-(MT(i,j)-Tin)/dx^2)-(-(MT(i+1,j)-MT(i,j))/dx^2)+...
            (-(MT(i,j)-MT(i,j-1))/dy^2)-(-(MT(i,j+1)-MT(i,j))/dy^2);
end
% C REGION
i=Nx;
for j=2:Ny-1
    MF(i,j)=(-(MT(i,j)-MT(i-1,j))/dx^2)-(-(Tout-MT(i,j))/dx^2)+...
            (-(MT(i,j)-MT(i,j-1))/dy^2)-(-(MT(i,j+1)-MT(i,j))/dy^2);
end
% D REGION
j=Ny;
for i=2:Nx-1
    MF(i,j)=(-(MT(i,j)-MT(i-1,j))/dx^2)-(-(MT(i+1,j)-MT(i,j))/dx^2)+...
            (-(MT(i,j)-MT(i,j-1))/dy^2)-(-Qtilde/dy);
end
% E REGION
j=1;
for i=2:Nx-1
     MF(i,j)=(-(MT(i,j)-MT(i-1,j))/dx^2)-(-(MT(i+1,j)-MT(i,j))/dx^2)+...
            0-(-(MT(i,j+1)-MT(i,j))/dy^2);
end
% F REGION
i=1;
j=Ny;
MF(i,j)=(-(MT(i,j)-Tin)/dx^2)-(-(MT(i+1,j)-MT(i,j))/dx^2)+...
            (-(MT(i,j)-MT(i,j-1))/dy^2)-(-Qtilde/dy);
%G REGION
i=Nx;
j=Ny;
MF(i,j)=(-(MT(i,j)-MT(i-1,j))/dx^2)-(-(Tout-MT(i,j))/dx^2)+...
            (-(MT(i,j)-MT(i,j-1))/dy^2)-(-Qtilde/dy);
% I REGION
i=1;
j=1;
MF(i,j)=(-(MT(i,j)-Tin)/dx^2)-(-(MT(i+1,j)-MT(i,j))/dx^2)+...
            0-(-(MT(i,j+1)-MT(i,j))/dy^2);
% H REGION
i=Nx;
j=1;
MF(i,j)=(-(MT(i,j)-MT(i-1,j))/dx^2)-(-(Tout-MT(i,j))/dx^2)+...
            0-(-(MT(i,j+1)-MT(i,j))/dy^2);

for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j;
        VF(k)=MF(i,j);
    end
end

VF=VF';

    
    
    
    
    