function z=CsignalDemo(P,w1,B,N,animate);
%Visual demo of  complex gaussian signals
% P is power, 
% w1 is angular center frequency
% N is number of samples
% B is rms bandwidth
% animate=1 gives temporal animation of signal, default=1
% 14.03.01  Hans Torp
% Aliasing frequency components included  23.03.04  H.T

if nargin<5,
    animate=1;
end;
if nargin<4,
    N=100;
end;

if nargin<1,%no input parameters specified
   P=1;
   w1=1;
   B=0.1;
   N=100;
   animate=1;
end; 


w=linspace(-pi,pi,N)';%angular frequency
t=1:N;%time in samples
m=t-N/2-1;%lag for autocorrelation function
m=-N+1:N-1;

%G=exp(-0.5*((w-w1)/B).^2);
G=exp(-0.5*((w-w1)/B).^2)+exp(-0.5*((w-w1-2*pi)/B).^2)+exp(-0.5*((w-w1+2*pi)/B).^2); %account for aliasing
G=P*G/mean(G);%Power spectrum
Zn=(randn(N,1)+i*randn(N,1))*sqrt(N/2);
Z=sqrt(fftshift(G)).*Zn;%shaping the dft of z
z=ifft(Z);%z is the sample process
if nargout>0,return;end;

R=fftshift(ifft(fftshift(sqrt(G))));R=xcorr(R);
figure(2);
subplot(3,1,1);plot(m,real(R),m,imag(R));grid;
title('Autocorrelation function');legend('real part','imaginary part');
subplot(3,1,2);plot(w,G);axis([-pi,pi,0,2*max(G)]);grid;
title('Power spectrum');
subplot(3,1,3);plot(w,10*log10(G));axis([-pi,pi,-40,20]);grid;
ylabel('Power spectrum [dB]');
xlabel('Angular frequency = w=2*pi*f/fs');
disp('Press space bar to continue');pause;
figure(1);
set(gcf,'Doublebuffer','on');%prevent flickering
subplot(2,2,2);plot(w,G);axis([-pi,pi,0,3*max(G)]);grid;
title('Power spectrum');
zm=2*sqrt(P);
if animate, n0=1, else n0=N;end;
for n=n0:N,
   subplot(2,2,1);
   plot(z(n),'*');
   axis('image');
   axis([-zm,zm,-zm,zm]);grid;
   xlabel('real part'); ylabel('imaginary part');title('Signal');
   subplot(2,2,3);
   t1=t(1:n);z1=z(1:n);
   plot(t1,real(z1),t1,imag(z1));
   axis([1,N,-zm,zm]);grid;legend('real part','imaginary part');
   pause(0.001);%to update display
end;
val=linspace(-zm,zm,10);
hx=hist(real(z),val);
hy=hist(imag(z),val);
subplot(2,2,4);plot(hx,val,hy,val);title('Histogram');
grid;legend('real part','imaginary part');

disp('Press space bar to continue');pause;
%power spectrum estimate
GN=fftshift(abs(fft(z)).^2)/N;
figure(2);
%GN=filter2(ones(50,1)/50,GN);
subplot(3,1,2);
plot(w,G,w,GN);axis([-pi,pi,0,3*max(G)]);grid;
legend('True spectrum','Estimated spectrum');
title('Power spectrum');
subplot(3,1,3);plot(w,10*log10(G),w,10*log10(GN));axis([-pi,pi,-40,20]);grid;
legend('True spectrum','Estimated spectrum');
ylabel('Power spectrum [dB]');
xlabel('Angular frequency = w=2*pi*f/fs');


