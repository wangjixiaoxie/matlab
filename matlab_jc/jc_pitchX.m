function [pitch_data]=jc_pitchX(shifted,N,OVERLAP,sigma,F_low,F_high,harms,filetype,returnAV)
%F_low and F_high are the boundary frequencies in which to look for a peak
%in the autocorrelation function (i.e. in which to calculate the pitch at
%each point in time).
%jc_pitchmat1024(shifted,1024,1010,2,2000,2700,[1 2 3],'obs0');

m=size(shifted);
NumberOfNotes=m(1);

%Parameters that one may want to change.
if strcmp(filetype,'obs0')
    SAMPLING=32000; 
else SAMPLING=44100;
end

%sonodvA program
t=-N/2+1:N/2;
sigma=(sigma/1000)*SAMPLING;
%Gaussian and first derivative as windows.
w=exp(-(t/sigma).^2);
timebins=floor((size(shifted,2)-(N))/(N-OVERLAP))+1;
freqbins=(N/2)+1;

Nyquist=SAMPLING/2;
step=Nyquist/freqbins; %This is the width (in ms) between bins.
mini=round(F_high/step);
maxi=round(F_low/step);
highest_harmonic=1; 

pitch_data=zeros(timebins,NumberOfNotes);
for found_note=1:NumberOfNotes
    
    sonos(found_note).data=abs(flipdim(spectrogram(shifted(found_note,:),w,OVERLAP,N,SAMPLING),1)); %gaussian windowed spectrogram

end

g=7;