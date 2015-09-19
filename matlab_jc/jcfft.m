function jcfft(x)
Fs = 1000; 
L = 1000; 
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(x,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2);
figure; plot(f,2*abs(Y(1:NFFT/2))) 