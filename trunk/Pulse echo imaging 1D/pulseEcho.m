%Pulse echo 1D imaging - simulation
% simulate received signal from an 1 D object defined by variation i acoustic impedance
%19.01.02  Hans Torp

%physical constants
c=1540;  %speed of sound

zmax=0.04;%max depth
dz=10e-6;% depth increment
z=0:dz:zmax;%depth axis
z=z';%make z into a vector

t=2/c*z;%time axis
dt=t(2)-t(1);
fs=1/dt;%sampling freq. 
tmax=t(end);

% calculate the total pulse echo response h = pel x pxd x hd x hxd ; where "x" means convolution
%define transmitted pulse
f0=6e6;
lambda=c/f0;
Tp=0.5e-6;%pulse length
tp=0:1/fs:Tp;tp=tp';
pel=sin(2*pi*f0*tp);% Electrical transmit pulse
figure(2);
plot(tp,pel);
%define transduce impulse response
fc=2.5e6;B=2.5e6;%center frequency and bandwidth of transducer
[bxd,axd]=butter(2,2*[fc-B/2,fc+B/2]/fs);%butterworth bandpass filter
txd=0:dt:3e-6;txd=txd';
impuls=zeros(size(txd));impuls(1)=1; %
hxd=filter(bxd,axd,impuls);
plot(txd,hxd);

pAc=conv(pel,hxd);%transmitted acoustical pulse
hd = [-1;1];% differensiation operator
h=conv(pAc,hd);%differensiation of Ac. impedance
h=conv(h,hxd);% transducer receiver response = transmit response hxd
plot(h);
%object defined by acoustic impedance Z as a function of depth
% Object 1: water/fat interface in depth z1
Zw=1.48; % [kg/m^2/s] acoustic impedance of water
Zf=1.37; % [kg/m^2/s] acoustic impedance of human fat tissue
z1=0.02;
Z=ones(size(z))*Zw;
Z(1+round(z1/dz):end)=Zf;
figure(1);subplot(3,1,1);
plot(z,Z);
%object defined by acoustic impedance Z as a function of depth
% Object 2: water/fat/water interface in depth z1, thickness th
Z=ones(size(z))*Zw;
z1=0.02;th=0.0005;
indLayer=round(z1/dz):round((z1+th)/dz);
Z(1+indLayer)=Zf;


s=conv(h,Z);%received signal
s=s(1:end-length(h)+1);
an=0.00025e-2;
sn=an*randn(size(s));%Gaussian white noise
s=s+sn;%add thermal noise

figure(1);
subplot(3,1,1);plot(z,Z);xlabel('depth z [m]');%axis('tight');
smax=max(abs(s));
subplot(3,1,2);plot(t,s);xlabel('time [s]');axis([t(1),t(end),-smax,smax]);
amp=abs(hilbert(s));%amplitude of echo signal s
%amp=abs((s));%amplitude of echo signal s
logamp=20*log10(amp);
subplot(3,1,3);plot(t,logamp);xlabel('time [s]');axis([t(1),t(end),-50,0]);


