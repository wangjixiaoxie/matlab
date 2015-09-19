function [x,psdvalues]=jcpsd2(pitch,fs)
% x is the period in milliseconds
% psdvalues is the POWER SPECTRAL DENSITY using matlab periodogram
normalized=pitch;%residuals(pitch);
for i=1:size(normalized,2)
    pitchN=normalized(:,i);
    [psdvalues(i,:),f]=periodogram(pitchN,[],[],fs);
end
    a=1000./f(2:end);
    x=[100;a]; % 100 is the DC component

