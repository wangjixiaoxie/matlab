function [CV,m]=CV730(pitchcurves)
for i=1:length(pitchcurves(1).pitches)
    for j=1:length(pitchcurves)
        slice(j)=pitchcurves(j).pitches(i);
    end
    m(i)=mean(slice);
    s=std(slice);
    CV(i)=s/m(i);
end

