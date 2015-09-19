function [medpitch,boomed]=jc_plotabunch721cbin(arrayfileD,f,xmax)

%notelength=850;
%Parameters
spacing=500;
note_length=6500;

sampling_rate=32000;

ymin=6500;
ymax=20000; 

xmin=0;
%xmax=length(arrayfileD(1).pitches);

%Smooth the data
for i=1:length(arrayfileD)
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

figure; hold on
for i=1:length(arrayfileD)
    h=xcov(avgnote(200:2200),smooth(i).smoothed(200:2200));
    [peak,index]=max(h);
    index=index+4500;
    shift(i)=(index-note_length)/4;
    %smooth
    processed=jc_pinterp(arrayfileD(i).pitches);
    medpitch(i)=median(processed(500:600));
    boomed(i)=std(processed(800:900));
    %shift them over to align
    k=1;
    for j=1:xmax
        align_shift=shift(i);
        shift_index=round(j-align_shift);
        if shift_index>0
            shifted(i,k)=processed(shift_index);
        else
            shifted(i,k)=ymin;
        end
        k=k+1;
    end
    plotshiftedD=shifted(i,:);
    ender=round(length(shifted)/23);
    for j=1:ender
        s(j)=std(plotshiftedD(j*20:j*20+60));
    end
    while j<40
        if s(j)<40 && s(j-1)< 40 && s(j-2)<40
            plotshiftedD=plotshiftedD(j*20:length(plotshiftedD));
            j=5000;
        else
        j=j+1;
        end
    end
    plotshiftedD=plotshiftedD+spacing*(i);
    plot(plotshiftedD); xlim([xmin xmax]); ylim([ymin ymax]); title('1-10')
end
