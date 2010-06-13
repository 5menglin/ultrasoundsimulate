
%Pulsed wave field from a focused line transducer
% Based on cwfieldfocusXd.m
%04.02.03  Hans Torp

f=5e6;%frekvens
c=1540;%lydhastighet
w=2*pi*f;%
k=w/c;%bølgetall
lambda=c/f%bølgelengde
%Transducer geometri
a=15e-3;%aperture
R=40e-3;%radius of curvature
dx0=lambda/2;%avstand mellom punktkilder på transducer
wAp=boxcar(round(a/dx0));%apodisation 
wAp=wAp(end/2:end);%cut away left half
%gating for pulsed wave
Tp=0.5e-6;%pulse length
tau=15e-6;%delay after pulse transmission
pw=1;%pw=1 turns range gating on

xmax=25e-3;% +/- størrelse i xretning
zmax=40e-3;%størrelse i zretning
Nx=345;%antall punkter i xretning
Nz=800;%antall punkter i zretning

xaxis=linspace(-xmax,xmax,Nx);
zaxis=linspace(0,zmax,Nz)';
x=ones(Nz,1)*xaxis;
z=zaxis*ones(1,Nx);

p0=0*x;%Akustisk felt ved tidspunkt t=0
figure(1);
set(gcf,'Doublebuffer','on');%prevent flickering
image(xaxis,zaxis,real(p0));
colormap(gray);axis('image');

p0=0*x;%setter akustisk felt til 0
n=1;
Np=length(p0(:));%total number of points
Ir=1:Np;
for x0=dx0/2:dx0:a/2,
   z0=R-sqrt(R^2-x0^2);
   r=sqrt((x-x0).^2+(z-z0).^2);%positiv x0
   if pw, Ir=find(abs(r/c-tau)<Tp/2);end;
   p0(Ir)=p0(Ir)+wAp(n)*exp(i*k*r(Ir))./r(Ir);
   r=sqrt((x+x0).^2+(z-z0).^2);%negativ x0
   if pw, Ir=find(abs(r/c-tau)<Tp/2);end;
   p0(Ir)=p0(Ir)+wAp(n)*exp(i*k*r(Ir))./r(Ir);
   n=n+1;
end;

image(xaxis,zaxis,0.05 *real(p0));axis image;

