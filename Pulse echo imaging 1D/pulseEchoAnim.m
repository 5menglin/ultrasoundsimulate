function pulsechoAnimate;
%Pulse echo 1D imaging - simulation
% simulate received signal from an 1 D object defined by variation in acoustic impedance
%29.01.05  Hans Torp

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
f0=2.5e6;
Tp=0.4e-6;%pulse length
tp=0:1/fs:Tp;tp=tp';
pel=sin(2*pi*f0*tp);% Electrical transmit pulse
figure(2);
Tpr=4e-6;%time-range 0 - Tpr for plot
subplot(3,2,1);plot(tp,pel);title('Electrical tx pulse');
xlim([0,Tpr]);
pause;
%define transducer impulse response
fc=2.5e6;B=2.2e6;%center frequency and bandwidth of transducer
[bxd,axd]=butter(2,2*[fc-B/2,fc+B/2]/fs);%butterworth bandpass filter
txd=0:dt:3e-6;txd=txd';
impuls=zeros(size(txd));impuls(1)=1; %
hxd=filter(bxd,axd,impuls);
subplot(3,2,2);plot(txd,hxd);title('transducer impulse response');
xlim([0,Tpr]);
pause;
 
pAc=conv(pel,hxd);%transmitted acoustical pulse
subplot(3,2,3);plot(pAc);title('transmitted acoustical pulse');
xlim([0,Tpr/dt]);
pause;
hd = 0.5*[-1;1];% differensiation operator to account for scattering
h=conv(hd,pAc);%differensiation of Ac. impedance
subplot(3,2,4);plot(h);title('reflected pulse');
xlim([0,Tpr/dt]);
pause;
h=conv(h,hxd);% transducer receiver response = transmit response hxd 
subplot(3,2,5);plot(h);title('received electrical signal');
xlim([0,Tpr/dt]);
pause;

%object defined by acoustic impedance Z as a function of depth
% Object 1: water/fat interface in depth z1
Zw=1.48; % [kg/m^2/s] acoustic impedance of water
Zf=1.37; % [kg/m^2/s] acoustic impedance of human fat tissue
z1=0.02;
Z=ones(size(z))*Zw;
Z(1+round(z1/dz):end)=Zf;

%object defined by acoustic impedance Z as a function of depth
%Object 2: water/fat/water interface in depth z1, thickness th
 Z=ones(size(z))*Zw;
 z1=0.02;th=0.01;
 indLayer=round(z1/dz):round((z1+th)/dz);
 Z(1+indLayer)=Zf;
 Z(1+indLayer)=Zf+0.01*randn(length(indLayer),1);

figure(1);clf;
subplot(3,1,1);
plot(z,Z);
if 1,
%animation of transmitted and scattered signal
figure(3);
ti=t/2;
pa=interp(pAc,2);
Np=length(pa);Tpac=Np*dt/2;
Nt=length(t);
Na=20;dta=Na*dt/2;
p0=zeros(round(Nt/Na),Nt);
pr=p0;
pr0=conv(h,Z/mean(Z));pr0=interp(pr0,2);
for k=1:size(p0,1),
    t1=k*dta;
    Iz=find( (t1-z/c)>0.5*dt & (t1-z/c)<Tpac );
    if ~isempty(Iz), 
        tp=t1-z/c;tp=tp(Iz);
        p0(k,Iz)=pAc(round(tp/dt));
    end;
    z1=c*t1;
    Iz=find( z<z1 );
    Ir1=2*k*Na;
    pr(k,Iz)=pr0(Ir1-length(Iz)+1:Ir1);
end;

for k=20:size(p0,1),
    figure(3);plot(z,Z,z,p0(k,:),z,pr(k,:));xlabel('Depth [m]');
    pause(0.05); 
    if k==20, pause;end
end;
pause;
%return;
end;


s=conv(h,Z/mean(Z));%received signal
s=s(1:end-length(h)+1);
an=0.05e-2;
sn=an*randn(size(s));%Gaussian white noise
s=s+sn;%add thermal noise

figure(1);
subplot(3,1,1);plot(z,Z);xlabel('depth z [m]');%axis('tight');
smax=max(abs(s));
subplot(3,1,2);plot(t,s);xlabel('time [s]');axis([t(1),t(end),-smax,smax]);
amp=abs(hilbert(s));%amplitude of echo signal s
logamp=20*log10(amp);
subplot(3,1,3);plot(t,logamp);xlabel('time [s]');axis([t(1),t(end),-80,0]);

%object defined by acoustic impedance Z as a function of depth
% Object 2: water/fat/water interface in depth z1, thickness th
% Z=ones(size(z))*Zw;
% z1=0.02;th=0.01;
% indLayer=round(z1/dz):round((z1+th)/dz);
% Z(1+indLayer)=Zf;
% 
