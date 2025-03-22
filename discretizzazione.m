clc
clear
delta=-0.5;
H0=2*abs(delta);
L=10*H0;
%fattore di finezza della mesh, numero intero positivo,
%all'aumentare del valore aumenta il numero di volumi finiti
fine=2;
Nx=(fine)*(H0*L)/(abs(delta));
dx=L/(Nx+1);
%per rispettare il rapporto tra i lati dei cateti il dy del triangolo è
%fissato
dy=dx*2*abs(delta)/L;
x=0:dx:L;
R=length(x)/2;
for i=1:R
H(i)=H0+delta*(2.*x(i))./L;
end
for i=(R+1):length(x)
H(i)=H0+delta-delta*((2.*x(i))./L-1);
end
Ny=-1+H./dy;
for i = 1:1:Nx
    for j = 1:1:Ny(i)
        k=sum(Ny(1:i-1))+j;
        a(k)=k;
    end
end
