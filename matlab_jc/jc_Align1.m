function shifted=jc_Align1(f)


note_length=6500;

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

for i=1:length(f)
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shifter(i)=(index-note_length);
end
    align_shift=shifter(i);
    beginning=1000-align_shift;
    shifted(i,:)=f(i).datt(beginning:beginning+4000);
end
