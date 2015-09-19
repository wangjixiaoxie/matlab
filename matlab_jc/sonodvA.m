function [sonogram,P]=sonodvA(s,SAMPLING,N,OVERLAP,sigma)

t=-N/2+1:N/2;
%Gaussian and first derivative as windows.

sigma=(sigma/1000)*SAMPLING;
w=exp(-(t/sigma).^2);
[sonogram,F,T,P]=spectrogram(s,w,OVERLAP,N,SAMPLING); %gaussian windowed spectrogram
sonogram=abs(flipdim(sonogram,1)); %this is the standard sonogram