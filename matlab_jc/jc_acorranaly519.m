function [covaried,avgcov]=jc_acorranaly519(arrayfile,f)


%Align
for i=1:length(arrayfile)
    [holder]=SmoothData(f(i).datt,32000,1);
    smooth(i).smoothed=holder;
end
for i=1:length(smooth(1).smoothed)
    sumnote(i)=0;
    for j=1:length(smooth)
        sumnote(i)=sumnote(i)+smooth(j).smoothed(i);
    end
    avgnote(i)=sumnote(i)/length(smooth);
end
for i=1:length(arrayfile)
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shift=(index-6400)/4;
    starting=round(280-shift);
    ending=starting+400;
    covaried(i).pitches=xcov(arrayfile(i).pitches(starting:ending));
end

%mean
for i=1:length(covaried(1).pitches)
    sumnote(i)=0;
    for j=1:length(covaried)
        sumnote(i)=sumnote(i)+covaried(j).pitches(i);
    end
    avgcov(i)=sumnote(i)/length(covaried);
end