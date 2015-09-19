function [Escape,notes]=jc_2contingA(pitchMat)
count=0;
for i=1:size(pitchMat,2)
    if mean(pitchMat(350:370,i))>mean(pitchMat(550:600,i))+80
        count=count+1;
        Escape(:,count)=pitchMat(:,i);
        notes(i)=10;
    end
end