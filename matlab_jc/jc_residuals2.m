function [residuals]=jc_residuals(pitch)
avgpitch=mean(pitch');
for i=1:size(pitch,1)
    residuals(i,:)=pitch(i,:)-avgpitch(i);
end