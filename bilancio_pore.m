function [FV]=bilancio_pore(t,CV)
parametri
alpha=alfa;
dx=1/(Nx-1);
dy=1/(Ny-1);

C=zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        k=(i-1)*Ny+j;
        C(i,j)=CV(k);
    end
end

% VECTOR FIELD FOR INNER-INNER NODES
for i=2:Nx-1
    for j=2:Ny-1
        aux=(C(i-1,j)-2*C(i,j)+C(i+1,j))/dx^2/alpha^2+(C(i,j-1)-2*C(i,j)+C(i,j+1))/dy^2;
        k=(i-1)*Ny+j;
        FV(k)=aux;
    end
end
% VECTOR FIELD FOR BOUNDARY NODES i=1 C(0,j)=1 
i=1;
for j=2:Ny-1
    aux=(1-2*C(i,j)+C(i+1,j))/dx^2/alpha^2+(C(i,j-1)-2*C(i,j)+C(i,j+1))/dy^2;
    k=(i-1)*Ny+j;
    FV(k)=aux;
end
% VECTOR FIELD FOR BOUNDARY NODES i=Nx C(Nx+1,j)=C(Nx,j)
i=Nx;
for j=2:Ny-1
    aux=(C(i-1,j)-2*C(i,j)+C(i,j))/dx^2/alpha^2+(C(i,j-1)-2*C(i,j)+C(i,j+1))/dy^2;
    k=(i-1)*Ny+j;
    FV(k)=aux;
end
% VECTOR FIELD FOR BOUNDARY NODES j=Ny+1 C(i,Ny+1)=C(i,Ny)/(1+Da*dy) 
j=Ny;
for i=2:Nx-1
    bv=C(i,j)/(1+Da*dy);
    aux=(C(i-1,j)-2*C(i,j)+C(i+1,j))/dx^2/alpha^2+(C(i,j-1)-2*C(i,j)+bv)/dy^2;
    k=(i-1)*Ny+j;
    FV(k)=aux;
end
% VECTOR FIELD FOR BOUNDARY NODES j=1 C(i,0)=C(i,1)/(1-Da*dy) 
j=1;
for i=2:Nx-1
    bv=C(i,j)/(1+Da*dy);
    aux=(C(i-1,j)-2*C(i,j)+C(i+1,j))/dx^2/alpha^2+(bv-2*C(i,j)+C(i,j+1))/dy^2;
    k=(i-1)*Ny+j;
    FV(k)=aux;
end
% VECTOR FIELD FOR BOUNDARY-CORNER NODES 
i=1;
j=1;
bv=C(i,j)/(1+Da*dy);
aux=(1-2*C(i,j)+C(i+1,j))/dx^2/alpha^2+(bv-2*C(i,j)+C(i,j+1))/dy^2;
k=(i-1)*Ny+j;
FV(k)=aux;
% ****************************
i=1;
j=Ny;
bv=C(i,j)/(1+Da*dy);
aux=(1-2*C(i,j)+C(i+1,j))/dx^2/alpha^2+(C(i,j-1)-2*C(i,j)+bv)/dy^2;
k=(i-1)*Ny+j;
FV(k)=aux;

i=Nx;
j=Ny;
bv=C(i,j)/(1+Da*dy);
aux=(C(i-1,j)-2*C(i,j)+C(i,j))/dx^2/alpha^2+(C(i,j-1)-2*C(i,j)+bv)/dy^2;
k=(i-1)*Ny+j;
FV(k)=aux;

i=Nx;
j=1;
bv=C(i,j)/(1+Da*dy);
aux=(C(i-1,j)-2*C(i,j)+C(i,j))/dx^2/alpha^2+(bv-2*C(i,j)+C(i,j+1))/dy^2;
k=(i-1)*Ny+j;
FV(k)=aux;


FV=FV';
end