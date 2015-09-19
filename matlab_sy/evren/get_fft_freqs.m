function ff=get_fft_freqs(nfft,fs);
% ff=get_fft_freqs(nfft,fs);
% NFFT is the # of points used in FFT
% fs is sampling freq in hz
% gives the negative freqs as well

len=nfft;
ff = [0:(nfft/2)]*fs/nfft;
ff = [ff,-ff(end-1:-1:2)];
return;
