function [Escape,notes]=jc_2conting(pitchMat)
count=0;
for i=1:size(pitchMat,2)
    if mean(pitchMat(345:370,i))>mean(pitchMat(545:570,i))+80
        count=count+1;
        Escape(:,count)=pitchMat(:,i);
        notes(i)=10;
    end
end