function yb = ExpBroaden(y,t)
% ExpBroaden(y,t) convolutes y by an exponential decay of time constant t
% by multiplying Fourier transforms and inverse transforming the result.
fy=fft(y);
a=exp(-[1:length(y)]./t);
fa=fft(a);
fy1=fy.*fa';           
yb=real(ifft(fy1))./sum(a);  
