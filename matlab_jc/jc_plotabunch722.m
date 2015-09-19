function [medpitch,boomed]=jc_plotabunch722(THarray,arrayfileD,f,xmax)

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
    processed3=jc_pinterp(THarray(i).pitches);

    %shift them over to align
    k=1;
    for j=1:xmax
        align_shift=shift(i);
        shift_index=round(j-align_shift);
        if shift_index>0
            shifted(i,k)=processed(shift_index);
            shifted3(i,k)=processed3(shift_index);
        else
            shifted(i,k)=ymin;
            shifted3(i,k)=ymin;
        end
        k=k+1;
    end
    plotshiftedD=shifted(i,:);
    plotshiftedD3=shifted3(i,:);
    ender=round(length(shifted)/23);
    for j=1:ender
        s(j)=std(plotshiftedD(j*20:j*20+60));
    end
    v=8;
    while v<40
        if s(v)<40 && s(v-1)< 40 && s(v-2)<40
            plotshiftedD=plotshiftedD(v*20-120:length(plotshiftedD));
            plotshiftedD3=plotshiftedD3(v*20-120:length(plotshiftedD3));
            v=5000;
        else
        v=v+1;
        end
    end
    medpitch(i)=mean(plotshiftedD3(200:350));
    plotshiftedD3=plotshiftedD3+spacing*(i);
    plot(plotshiftedD3); xlim([xmin xmax]); ylim([ymin ymax]); title('1-10')
end
