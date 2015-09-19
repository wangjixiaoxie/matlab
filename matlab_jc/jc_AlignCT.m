function [shifted]=jc_AlignCT(f,toff,ind)


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
        if isempty(find(ind==i))

            beginning=round(toff(count)*32-3200);
            %if abs(align_shift(i))<600
            if count==26
                beginning=500;
            end
                shifted(count,:)=zeros(1,10500);
                %shifted(count,2200-beginning:2200-beginning+6599)=f(i).datt(1:6600);
                %count=count+1;
                begger(i)=beginning;
            %else 
             %   i
            %end
            count=count+1;
        end
        
end
g=8;




