function [CV]=CV731window(pitchcurves)
b=[4,8,16,32,64,128];
for i=1:length(b)
    window=b(i);
    for j=1:length(pitchcurves)
        slice=pitchcurves(j).pitches(200-window:200+window);
        mu(j)=mean(slice);
    end
    m=mean(mu);
    s=std(mu);
    CV(i)=s/m;
end

