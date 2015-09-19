function shifted=jc_AlignPart(f)


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
    h=xcov(avgnote(1:5000),smooth(i).smoothed(1:5000));
    [peak,index]=max(h);
    shifter=(index-length(avgnote(1:5000)));
    
    align_shift(i)=shifter;
    beginning=800-align_shift(i);
    shifted(i,:)=f(i).datt(beginning:beginning+12000);
end


