function [averaged,SE,change]=jc_getavstd523(arrayfile,f)

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
    processed(i).pitches=arrayfile(i).pitches(starting:ending);
    left=arrayfile(i).pitches(starting);
    right=arrayfile(i).pitches(ending);
    first(i)=left;
    change(i)=(right-left);
end

%Take the average and standard deviation of the aligned data
for i=1:length(processed(1).pitches)
    for j=1:length(processed)
        matr(i,j)=processed(j).pitches(i);
    end
    averaged(i)=mean(matr(i,:));
    standev(i)=std(matr(i,:));
    SE(i)=standev(i)/sqrt(length(processed));
end

