clc
clear
%parametro geometrico variabile
delta=0;
%parametri geometrici fissi del canale
H0=1;
L=10;
%numero intero da scegliere
Nx0=500;
%correzione al numero di volumi per garantire di avere Ny intero
Nx=Nx0*(abs(delta))/H0;
if delta == 0
    Nx = Nx0;
end
%definizione della tasselazione in direzione trasversale
for i=1:Nx
    Ny(i)=Nx0+i;
end
for i=1:Nx
    Ny(2*Nx+1-i)=Nx0+i;
end
Nx=2*Nx;
%definizione delle dimensioni dei volumi finiti
dx=L/Nx;
dy=abs(delta)*(L/2)*dx;
if delta == 0
    dy = dx;
end
%definizione del sistema di coordinate
x=0:dx:L;
maxLength = length(0:dy:max(H0, H0+delta));
y = zeros(maxLength, Nx); % Preallocate matrix
for i = 1:(Nx/2)
    temp = 0:dy:(H0 + delta * 2 * i / Nx);
    y(1:length(temp), i) = temp; % Assign correctly
end
for i = 1:(Nx/2)
    temp = 0:dy:(H0 + delta * 2 * i / Nx);
    y(1:length(temp), Nx+1-i) = temp; % Assign correctly
end
%passaggio da un array a due indici a un vettore
for i = 1:1:Nx
    for j = 1:1:Ny(i)
        k=sum(Ny(1:i-1))+j;
        a(k)=k;
    end
 end
