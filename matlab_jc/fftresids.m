function fftresids(y)
a=size(y);
for i=1:a(2)
m=y(:,i);
L=length(m);
Fs=30;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(m,NFFT)/L;
Z(i,:)=2*abs(Y(1:NFFT/2));
end
for j=1:length(Z)
    K(j)=mean(Z(:,j));
end
% Plot single-sided amplitude spectrum.
plot(K) 