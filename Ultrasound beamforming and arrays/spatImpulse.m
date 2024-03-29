%spatial impulse response from a discrete aperture with points (xa,za) in point x,z
%define aperture
a=16e-3;
R=50e-3;
da=0.02e-3;
xa=-a/2:da:a/2;
za=R-sqrt(R^2-xa.^2);
x=(-25:0.25:25)*1e-3;%x-positions
z=50e-3;

c=1540;
fs=40e6;%temp. sampling frequency
T=40e-6;%time axis 0 - T

dt=1/fs;
%t=0:dt:T;t=t';
%t=0.75*T:dt:T;t=t';
t=0.75*T:dt:T;t=t';
Na=length(xa);
h=zeros(length(t),length(x));
for m=1:length(x),
    h1=0*t;
    for n=1:Na,
        r=sqrt((x(m)-xa(n))^2+(z-za(n))^2);
        tau=r/c;
        ind=1+round((tau-t(1))/dt);
        ind=ind-1:ind+1;
        h1(ind)=h1(ind)+1/r;
    end;
    %figure(1);plot(t,h1);pause(0.05);
    h(:,m)=h1;
end;
h=h/max(h(:));  %规一化的操作。
figure(2);subplot(2,1,1);mesh(x,t,h);
pause;%return;
%transmitted pulse
f0=3e6;
Tp=0.33e-6;%pulse length
tp=0:1/fs:Tp;tp=tp';
pel=sin(2*pi*f0*tp);% Electrical transmit pulse
figure(1);plot(tp,pel)
%define transduce impulse response
fc=3e6;B=2e6;%center frequency and bandwidth of transducer
[bxd,axd]=butter(2,2*[fc-B/2,fc+B/2]/fs);%butterworth bandpass filter
txd=0:dt:3e-6;txd=txd';
impuls=zeros(size(txd));impuls(1)=1; %
hxd=filter(bxd,axd,impuls);
pAc=conv(pel,hxd);%transmitted acoustical pulse
plot(pAc)
hp=filter(pAc,1,h);
hp=hp/max(hp(:));%normaliserer
for m=1:length(x)/2,
    figure(1);plot(t,h(:,m),t,hp(:,m));
    axis([3e-5,4e-5,-1,1]);
    pause;
end;
figure(2);subplot(2,1,2);pcolor(x,t,20*log10(abs(hp)));
caxis([-40,0]);shading interp;axis ij;colorbar;

%spatial frequency response
H=fft(h);H=H/max(H(:));
Hp=fft(hp);Hp=Hp/max(Hp(:));
Nf=length(t);
f=(0:(Nf-1))/Nf*fs;
pause;
figure(3);
pcolor(x,f,20*log10(abs(Hp)));
caxis([-50,0]);shading interp;axis ij;colorbar;
axis([x(1),x(end),0,10e6]);
return;
%aperturefunction in k-space
H=fft2(hp);H=H/max(H(:));
H=abs(fftshift(H));
%H=min(0.01,H);
%mesh(H);

figure(3);
imagesc(H);
pcolor(abs(H));
caxis([0,1]);shading interp;axis ij;
%axis([x(1),x(end),0,10e6]);


