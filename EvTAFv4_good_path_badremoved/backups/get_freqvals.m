function FreqVals=get_freqvals(NFFT,FSample);
%FreqVals=get_freqvals(NFFT,FSample);
% this is a function i wrote in LabVIEW to get the frequency values for each freq bin of the
% FFT done in LV hopefully it's right and matches what labview gets

NN=fix(round(NFFT/2));
DF=FSample./2./NN;
FreqVals=[0:NN-1]*DF;
return;
