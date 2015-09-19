function [x,psdmedian]=jc_psd(pitch)
%%%% Takes as input a clean region of pitch curves e.g.jc_psd(pitchUDpre(450:800,:))

%%%% Take the residuals
residuals=jc_residuals(pitch);

%%%% Take the autocorrelation function of each residual
autocorr=jc_xcorr(residuals);

%%%% Take the PSD (FFT of each autocorrelation function)
%Fs=8; % per millisecond --- 8 for STFT, x for hilbert
for i=1:size(autocorr,2)
    acorrinstance=autocorr(:,i);
    L=length(acorrinstance);
    NFFT=2^nextpow2(L);
    % no normalization for bandwidth
    psda=fft(acorrinstance,NFFT);
    psdvalues(:,i)=2*abs(psda(1:NFFT/2));
end

%%%% Take the mean of the PSDs -- an issue that needs resolution is why I
%%%% get such crappy results from taking the median.
psdmedian=median(psdvalues');

% x-axis calibration -- period of sine waves
% http://www.dsprelated.com/showmessage/47831/1.php
fs=32000;
long=length(psdmedian);
x(1)=100; % DC component
for i=1:long-1
    x(i+1)=(1/(fs*(i/long)))*1000;
end