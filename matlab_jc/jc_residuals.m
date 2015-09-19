function [residuals]=jc_residuals(pitch)
avgpitch=mean(pitch');
for i=1:size(pitch,1)  % this has been a correction
    residuals(i,:)=pitch(i,:)./avgpitch(i)-1;
end