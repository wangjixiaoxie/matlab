function [means]=jc_getmeans526(arrayfile,f)

%Smooth the data
for i=1:length(arrayfile)
    [holder]=SmoothData(f(i).datt,32000,1);
    smooth(i).smoothed=holder;
end

%Get average smoothed note
for i=1:length(smooth(1).smoothed)
    sumnote(i)=0;
    for j=1:length(smooth)
        sumnote(i)=sumnote(i)+smooth(j).smoothed(i);
    end
    avgnote(i)=sumnote(i)/length(smooth);
end

%Align the data
for i=1:length(arrayfile)
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shift=(index-6400)/4;
    starting=round(280-shift);
    ending=starting+400;
    for n=1:3
        processed(i).freqbinest(n,:)=arrayfile(i).freqbinest(n,starting:ending);
    end
end

for i=1:length(processed)
    for n=1:3
        means(i,n)=mean(processed(i).freqbinest(n));
    end
end
