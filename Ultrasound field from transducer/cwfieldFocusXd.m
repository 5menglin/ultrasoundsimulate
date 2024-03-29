
% Continuous wave from fokused line-transducer
%24.01.01  Hans Torp
%Modified to include apodisation 24.01.01  Hans Torp

f=2e6;%frekvens 频率
c=1540;%lydhastighet 声速
w=2*pi*f;%角速率
k=w/c;%b鴏getall
lambda=c/f%b鴏gelengde 人 波长
animate=1;%是否进行动画显示

xmax=25e-3;% +/- st鴕relse i xretning
zmax=50e-3;%st鴕relse i zretning
Nx=500;%antall punkter i xretning
Nz=500;%antall punkter i zretning

xaxis=linspace(-xmax,xmax,Nx);
zaxis=linspace(0,zmax,Nz)';
x=ones(Nz,1)*xaxis;             % 对于每一列,X的位置
z=zaxis*ones(1,Nx);             % 对于每一行，Y的位置

p0=0*x;%Akustisk felt ved tidspunkt t=0
figure(1);
set(gcf,'Doublebuffer','on');%prevent flickering
image(xaxis,zaxis,real(p0));
colormap(gray);axis('image');
dt=1/f/10;T=1/f;%for animation

%Transducer geometri
a= 20e-3;%aperture
R= 35e-3;%radius of curvature
dx0=2*lambda%avstand mellom punktkilder p?transducer
wAp=boxcar(round(a/dx0));%apodisation %矩形窗
%wAp=hamming(round(a/dx0));%apodisation 
wAp=wAp(end/2:end);%cut away left half


p0=0*x;%setter akustisk felt til 0
n=1;
for x0=dx0/2:dx0:a/2,
   z0=R-sqrt(R^2-x0^2);
   r=sqrt((x-x0).^2+(z-z0).^2);%positiv x0
   p0=p0+wAp(n)*exp(i*k*r)./r;
   r=sqrt((x+x0).^2+(z-z0).^2);%negativ x0
   p0=p0+wAp(n)*exp(i*k*r)./r;
   n=n+1;
   crossover=(2*x0)^2/lambda;
   if animate,%animation of wave propagation
       T=1/f;
       for t=0:dt:T,
           p=p0*exp(-i*w*t);
           image(xaxis,zaxis,0.04*real(p));
           axis('image');pause(0.001);
       end;
   end;
   %pause
end;
image(xaxis,zaxis,0.2*real(p0));
axis('image');
return;
if animate,%animation of wave propagation
    dt=1/f/10;T=4/f;
    for t=0:dt:T,
        p=p0*exp(-i*w*t);
        image(xaxis,zaxis,0.04*real(p));
        axis('image');pause(0.001);
    end;
end;

%str錶eprofil
% Pr鴙 ?sammenlign med de analytiske uttrykkene i N. Wright's kompendium!
z0=1;
pmax=max(max(abs(p0(100:end,:))));
while z0>0,
    disp('Str錶eprofil: Klikk i bildet');
    figure(1);
    [x0,z0]=ginput(1);
    if z0>0,
        dz=zaxis(2)-zaxis(1);
        z0ind=1+round(z0/dz);
        dBprofil=20*log10(abs(p0(z0ind,:)/pmax));
        figure(2);plot(xaxis,dBprofil);axis([-xmax,xmax,-50,0]);
        grid;
    end;
end;

