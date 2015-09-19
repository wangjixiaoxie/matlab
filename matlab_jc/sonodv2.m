function [sonogram]=sonodv(s,SAMPLING,N,OVERLAP,sigma)

sonogram=spectrogram(s,64,60,64,32000)+eps; %gaussian windowed spectrogram
sonogram=abs(flipdim(sonogram,1)); %this is the standard sonogram