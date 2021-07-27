%keplerian turbulence synthetic model (classic)

%wave number parameters
N=512;
kmin=5;
kmax=50;

%cell parameters
M=64;
Lx=2;
Ly=2;
Lz=2;

for j=1:M
    x(j)=-Lx/2+(j-1)*Lx/(M-1);
    y(j)=-Ly/2+(j-1)*Ly/(M-1);
    z(j)=-Lz/2+(j-1)*Lz/(M-1);
end

%wave vector amplitudes
kamp(1)=kmin;
kamp(N)=kmax;
for j=2:(N-1)
    kamp(j)=kamp(1)*(kamp(N)/kamp(1))^((j-1)/(N-1));
end

%energy dissipation
eps=(1.5/((kmin)^(-2/3)-(kmax)^(-2/3)))^(3/2);

%Kolmogorov energy spectrum and frequencies
for j=1:N
    E(j)=1.5*(eps)^(2/3)*kamp(j)^(-5/3);
    w(j)=0.5*sqrt(kamp(j)^3*E(j));
end

%amplitudes of vectors a and b
deltak(1)=(kamp(2)-kamp(1))/2;
deltak(N)=(kamp(N)-kamp(N-1))/2;
for i=2:(N-1)
    deltak(i)=(kamp(i+1)-kamp(i-1))/2;
end

for i=1:N
    aamp(i)=sqrt(2*E(i)*deltak(i));
    bamp(i)=sqrt(2*E(i)*deltak(i));
end

%direction of vectors k,a,b
for i=1:N
    %k
fi=2*pi*rand(1);
tetta=pi-pi*rand(1);
kedx(i)=sin(tetta)*sin(fi);
kedy(i)=sin(tetta)*cos(fi);
kedz(i)=cos(tetta);
kx(i)=kedx(i)*kamp(i);
ky(i)=kedy(i)*kamp(i);
kz(i)=kedz(i)*kamp(i);
    %a
aedx1(i)=1-2*rand(1);
aedy1(i)=1-2*rand(1);
aedz1(i)=-(kedx(i)*aedx1(i)+kedy(i)*aedy1(i))/kedz(i);
aedamp1(i)=sqrt(aedx1(i)^2+aedy1(i)^2+aedz1(i)^2);
aedx(i)=aedx1(i)/aedamp1(i);
aedy(i)=aedy1(i)/aedamp1(i);
aedz(i)=aedz1(i)/aedamp1(i);
ax(i)=aedx(i)*aamp(i);
ay(i)=aedy(i)*aamp(i);
az(i)=aedz(i)*aamp(i);
    %b
bedx1(i)=1-2*rand(1);
bedy1(i)=1-2*rand(1);
bedz1(i)=-(kedx(i)*bedx1(i)+kedy(i)*bedy1(i))/kedz(i);
bedamp1(i)=sqrt(bedx1(i)^2+bedy1(i)^2+bedz1(i)^2);
bedx(i)=bedx1(i)/bedamp1(i);
bedy(i)=bedy1(i)/bedamp1(i);
bedz(i)=bedz1(i)/bedamp1(i);
bx(i)=bedx(i)*bamp(i);
by(i)=bedy(i)*bamp(i);
bz(i)=bedz(i)*bamp(i);
end

clear kedy kedx kedz deltak aamp bamp aedx1 aedy1 aedz1 bedx1 bedy1 bedz1 aedx aedy aedz bedx bedy bedz aedamp1 bedamp1;
%synthetic velocity field in the real space
dt=0.01;
im(1,1,1,501) = 0; 
hF=figure; 

for j=1:M
    for l=1:M
        for m=1:M
            for d=1:N                
vx1(j,l,m,d)=(ax(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+bx(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
ddxvx1(j,l,m,d)=(-ax(d)*kx(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+bx(d)*kx(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
ddyvx1(j,l,m,d)=(-ax(d)*ky(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+bx(d)*ky(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
ddzvx1(j,l,m,d)=(-ax(d)*kz(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+bx(d)*kz(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
            end
        end
    end
    j
end

vx=sum(vx1,4);
ddxvx=sum(ddxvx1,4);
ddyvx=sum(ddyvx1,4);
ddzvx=sum(ddzvx1,4);
clear vx1 ddxvx1 ddzvx1;

for j=1:M
    for l=1:M
        for m=1:M
            for d=1:N                
vy1(j,l,m,d)=(ay(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+by(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
ddxvy1(j,l,m,d)=(-ay(d)*kx(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+by(d)*kx(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
ddyvy1(j,l,m,d)=(-ay(d)*ky(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+by(d)*ky(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
ddzvy1(j,l,m,d)=(-ay(d)*kz(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+by(d)*kz(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
            end
        end
    end
    j
end

vy=sum(vy1,4);
ddxvy=sum(ddxvy1,4);
ddyvy=sum(ddyvy1,4);
ddzvy=sum(ddzvy1,4);
clear vy1 ddxvy1 ddzvy1;

for j=1:M
    for l=1:M
        for m=1:M
            for d=1:N                
vz1(j,l,m,d)=(az(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+bz(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
ddxvz1(j,l,m,d)=(-az(d)*kx(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+bz(d)*kx(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
ddyvz1(j,l,m,d)=(-az(d)*ky(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+bz(d)*ky(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
ddzvz1(j,l,m,d)=(-az(d)*kz(d)*sin(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m))+bz(d)*kz(d)*cos(kx(d)*x(j)+ky(d)*y(l)+kz(d)*z(m)));
            end
        end
    end
    j
end

vz=sum(vz1,4);
ddxvz=sum(ddxvz1,4);
ddyvz=sum(ddyvz1,4);
ddzvz=sum(ddzvz1,4);
clear vy1 ddxvy1 ddzvy1;

for j=1:M
    for l=1:M
        for m=1:M
Density(j,l,m)=0.02*(ddxvx(j,l,m)^2+ddyvy(j,l,m)^2+ddzvz(j,l,m)^2+2*(ddxvy(j,l,m)*ddyvx(j,l,m)+ddxvz(j,l,m)*ddzvx(j,l,m)+ddyvz(j,l,m)*ddzvy(j,l,m))-3*ddyvx(j,l,m));
        end
    end
end

%Density at the midplane z=0:
for i=1:M
    for j=1:M
        Density0(i,j)=Density(i,j,0);
    end
end

%vyrez=0.25*sum(vy1,4)+u0;

%kinetic energy distribution at the physical space
for j=1:M
    for l=1:M
        for m=1:M
            EE(j,l,m)=(1/2)*(vx(j,l,m)^2+vy(j,l,m)^2+vz(j,l,m)^2);
        end
    end
end

isocaps(x,y,z,vx,-10)
view(20,20);
title(['Ux(x,y,z,0)']);
colorbar('vert');
caxis([-6, 6]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;


hF=figure; 
isocaps(x,y,z,vy,-10)
view(20,20);
title(['Uy(x,y,z,0)']);
colorbar('vert');
caxis([-6, 6]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;

hF=figure; 
isocaps(x,y,z,vz,-10)
view(20,20);
title(['Uz(x,y,z,0)']);
colorbar('vert');
caxis([-6, 6]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;

hF=figure; 
isocaps(x,y,z,Density,-200)
view(20,20);
title(['Density(x,y,z)']);
colorbar('vert');
caxis([-30, 30]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;


for j=1:M
    z0(j)=0;
end
 
hF=figure; 
isocaps(x,y,z0,Density,-200)
view(0,90);
title(['Density(x,y,z)']);
colorbar('vert');
caxis([-30, 30]);
xlabel('x'); 
ylabel('y');
zlabel('z'); 
grid on;


