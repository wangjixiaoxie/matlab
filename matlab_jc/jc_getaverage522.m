function [avgup,avgdown,avgoverall,j,k,change]=jc_getaverage522(arrayfile,f)

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
    fifty(i)=arrayfile(i).pitches(starting+350);
    change(i)=(right-left);
end
avg=mean(fifty);
j=1;
k=1;
for i=1:length(arrayfile)
    if fifty(i)>avg+1*std(fifty)
        avup(j).pitches=processed(i).pitches;
        j=j+1;
    else
        if fifty(i)<avg-1*std(fifty)
        avdown(k).pitches=processed(i).pitches;
        k=k+1;
        end
    end
end
        
%Take the average of the aligned data
for i=1:length(avup(1).pitches)
    summed(i)=0;
    for j=1:length(avup)
        summed(i)=summed(i)+avup(j).pitches(i);
    end
    avgup(i)=summed(i)/length(avup);
end
for i=1:length(avdown(1).pitches)
    summed(i)=0;
    for j=1:length(avdown)
        summed(i)=summed(i)+avdown(j).pitches(i);
    end
    avgdown(i)=summed(i)/length(avdown);
end
for i=1:length(processed(1).pitches)
    summed(i)=0;
    for j=1:length(processed)
        summed(i)=summed(i)+processed(j).pitches(i);
    end
    avgoverall(i)=summed(i)/length(processed);
end

