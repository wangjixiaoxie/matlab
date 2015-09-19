function [align_shift]=jc_align_shift(f)




sampling_rate=32000;


%Smooth the data
for i=1:length(f)
    [holder]=SmoothData(f(i).datt,sampling_rate,1);
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
count=1;
vals=getvals(f,1,'TRIG');
for i=1:length(f)
    if vals(i,3)~=0
        h=xcov(avgnote(3000:7500),smooth(i).smoothed(3000:7500));
        [peak,index]=max(h);
        shifter=(index-length(avgnote(3000:7500)));
        align_shift(count)=shifter;
        count=count+1;
    end
end
g=8;




