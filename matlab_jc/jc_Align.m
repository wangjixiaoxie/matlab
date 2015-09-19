function [align_shift]=jc_Align(f)


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
count=1;
for i=1:length(f)
    if i~=8 && i~=50 && i~=63 && i~=69
    
        h=xcov(avgnote(3000:7500),smooth(i).smoothed(3000:7500));
        [peak,index]=max(h);
        shifter=(index-length(avgnote(3000:7500)));

        align_shift(count)=shifter;
        %beginning=600-align_shift(i);
        %if abs(align_shift(i))<600
            %shifted(count,:)=f(i).datt(beginning:beginning+7000);
            %count=count+1;
        %else 
         %   i
        %end
        count=count+1;
    end
end
g=8;




