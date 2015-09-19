function h=jc_fft(g)
ll=length(g);
nfft=2^nextpow2(ll);
Y=fft(g,nfft);
f=8000/2*linspace(0,1,nfft/2);
h=2*abs(Y(1:nfft/2));
