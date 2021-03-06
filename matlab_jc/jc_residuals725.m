function [xcresids,normalized,avgnote]=jc_getresiduals616(arrayfile,f,xmax,starting,ending)
%800,280,680 or 1000,400,800
note_length=6400;
sampling=44100;



ymin=0;




%Smooth the data
for i=1:length(arrayfile)
    [holder]=SmoothData(f(i).datt,sampling,1);
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
    file=arrayfile(i).pitches;
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shift(i)=(index-note_length)/4;
    %shift them over to align
    k=1;
    for j=1:xmax
        align_shift=shift(i);
        shift_index=round(j-align_shift);
        if shift_index>0
            shifted(i,k)=file(shift_index);
        else
            shifted(i,k)=ymin;
        end
        k=k+1;
    end
    processed(i).pitches=shifted(i,starting:ending);
end



%Get the residuals
for i=1:length(processed(1).pitches)
    for j=1:length(processed)
        matr(i,j)=processed(j).pitches(i);
    end
    averaged(i)=mean(matr(i,:));
    for j=1:length(processed)
        normalized(i,j)=matr(i,j)/averaged(i);
        normalized(i,j)=normalized(i,j)-1;
    end
end

%Take the cross-correlation of the residuals for each trace.
for j=1:length(processed)
    xcr(:,j)=xcov(normalized(:,j));
end


for i=1:length(xcr)
    xcresids(i)=mean(xcr(i,:));
end