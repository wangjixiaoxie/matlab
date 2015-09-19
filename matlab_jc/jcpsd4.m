function [x,psdvalues]=jcpsd3(pitch,fs)
% x is the period in milliseconds
% psdvalues is the POWER SPECTRAL DENSITY using matlab periodogram
pitchN=pitch-mean(pitch);
    [psdvalues,f]=periodogram(pitchN,[],[],fs);
    a=1000./f(2:end);
    x=[100;a]; % 100 is the DC component

