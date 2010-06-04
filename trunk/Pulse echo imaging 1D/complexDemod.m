%Complex demodulation- statistical properties
% 18.01.07 Hans Torp
T=2500e-6;
fs=50e6;
f0=2e6;
B=0.5e6;


dt=1/fs;
t=0:1/fs:T;
s=randn(size(t));%Gaussian random noise
%create a bandpass signal 
[b,a]=butter(2,2*[f0-B/2,f0+B/2]/fs);%butterworth bandpass filter
s=filter(b,a,s);
figure(1);
subplot(2,2,1);plot(t,s);
subplot(2,2,2);hist(s,30);

%complex demodulation
sp=hilbert(s);
z=sp.*exp(-i*2*pi*f0*t);
subplot(2,2,3);plot(t,s,t,real(z),t,imag(z),t,abs(z));zoom xon;
legend('s','real','imag','abs(z)');
subplot(2,2,4);hist([real(z);imag(z)]',30);
legend('real','imag');
figure(2);plot(z,'.');axis('equal');

P=mean(abs(z).^2);
P=mean(z.*conj(z))
Rzz=mean(z.*z)
abs(Rzz)/P

abs(mean(z.*z))/mean(z.*conj(z))


