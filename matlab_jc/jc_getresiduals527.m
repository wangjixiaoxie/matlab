function [residuals,residPLOT]=jc_getresiduals527(arrayfile,f)

note_length=4096;
sampling=44100;


xmax=500;
ymin=0;
starting=50;
ending=250;



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
    residPLOT(i)=std(normalized(i,:));
end

%Put the residuals in a 1xn matrix
k=0;
for i=1:length(processed(1).pitches)
    for j=1:length(processed)
        k=k+1;
        residuals(k)=normalized(i,j);
    end
end