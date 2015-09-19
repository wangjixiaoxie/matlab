function [x,psdvalues]=jcpsd(pitchN,fs)
% x is the period in milliseconds
% psdvalues is the POWER SPECTRAL DENSITY using matlab periodogram

[psdvalues,f]=periodogram(pitchN,[],[],fs);
a=1000./f(2:end);
x=[100;a]; % 100 is the DC component