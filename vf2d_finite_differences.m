function [B]=vf2d(t,ck)
parameters2d_finite_differences

for i=1:Nx
    for j=1:Nz
        k=(i-1)*Nz+j;
        c(i,j)=ck(k);
    end
end

dx=(L/H)/(Nx+1);
dz=1/(Nz+1);
for j=1:Nz
    z(j)=j*dz;
    f(j)=6*z(j)*(1-z(j));
end

% internal internal points
for i=2:Nx-1
    for j=2:Nz-1
     vf(i,j)=(c(i+1,j)-2*c(i,j)+c(i-1,j))/Pe/dx^2+(c(i,j+1)-2*c(i,j)+c(i,j-1))/Pe/dz^2-...
         f(j)*(c(i,j)-c(i-1,j))/dx-Da*c(i,j)^n;
    end
end
% boundary points A
j=1;
for i=2:Nx-1
    vf(i,j)=(c(i+1,j)-2*c(i,j)+c(i-1,j))/Pe/dx^2+(c(i,j+1)-2*c(i,j)+c(i,j))/Pe/dz^2-...
         f(j)*(c(i,j)-c(i-1,j))/dx-Da*c(i,j)^n;
end
% boundary points C
j=Nz;
for i=2:Nx-1
    vf(i,j)=(c(i+1,j)-2*c(i,j)+c(i-1,j))/Pe/dx^2+(c(i,j)-2*c(i,j)+c(i,j-1))/Pe/dz^2-...
         f(j)*(c(i,j)-c(i-1,j))/dx-Da*c(i,j)^n;
end
% boundary points B
i=Nx;
for j=2:Nz-1
    vf(i,j)=(c(i,j)-2*c(i,j)+c(i-1,j))/Pe/dx^2+(c(i,j+1)-2*c(i,j)+c(i,j-1))/Pe/dz^2-...
         f(j)*(c(i,j)-c(i-1,j))/dx-Da*c(i,j)^n;
end
i=1;
% boundary points D
for j=2:Nz-1
    vf(i,j)=(c(i+1,j)-2*c(i,j)+c0)/Pe/dx^2+(c(i,j+1)-2*c(i,j)+c(i,j-1))/Pe/dz^2-...
         f(j)*(c(i,j)-c0)-Da*c(i,j)^n;
end
% boundary point E
i=1;
j=1;
vf(i,j)=(c(i+1,j)-2*c(i,j)+c0)/Pe/dx^2+(c(i,j+1)-2*c(i,j)+c(i,j))/Pe/dz^2-...
      f(j)*(c(i,j)-c0)/dx-Da*c(i,j)^n;
% boundary point F
i=Nx;
j=1;
vf(i,j)=(c(i,j)-2*c(i,j)+c(i-1,j))/Pe/dx^2+(c(i,j+1)-2*c(i,j)+c(i,j))/Pe/dz^2-...
         f(j)*(c(i,j)-c(i-1,j))/dx-Da*c(i,j)^n;
% boundary point G
i=Nx;
j=Nz;
vf(i,j)=(c(i,j)-2*c(i,j)+c(i-1,j))/Pe/dx^2+(c(i,j)-2*c(i,j)+c(i,j-1))/Pe/dz^2-...
         f(j)*(c(i,j)-c(i-1,j))/dx-Da*c(i,j)^n;     
% boundary point H
i=1;
j=Nz;
vf(i,j)=(c(i+1,j)-2*c(i,j)+c0)/Pe/dx^2+(c(i,j)-2*c(i,j)+c(i,j-1))/Pe/dz^2-...
         f(j)*(c(i,j)-c0)/dx-Da*c(i,j)^n;

for i=1:Nx
    for j=1:Nz
        k=(i-1)*Nz+j;
        B(k)=vf(i,j);
    end
end     
B=B';
     

     
     
     




       

    