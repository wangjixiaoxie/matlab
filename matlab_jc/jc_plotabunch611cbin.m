function shift=jc_plotabunch611cbin(arrayfileD,f,xmax)

%notelength=850;
%Parameters
spacing=180;
note_length=6400;

sampling_rate=32000;

ymin=2400;
ymax=4400; 

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
subplot(131); hold on
for i=1:10
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shift(i)=(index-note_length)/4;
    %smooth
    processed=jc_pinterp(arrayfileD(i).pitches);
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
    plotshiftedD=shifted(i,:)+spacing*(i);
    plot(plotshiftedD); xlim([xmin xmax]); ylim([ymin ymax]); title('1-10')
end

subplot(132); hold on
for i=11:20
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shift(i)=(index-note_length)/4;
    %smooth
    processed=jc_pinterp(arrayfileD(i).pitches);
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
    plotshiftedD=shifted(i,:)+spacing*(i-10);
    plot(plotshiftedD); xlim([xmin xmax]); ylim([ymin ymax]); title('10-20')
end

subplot(133); hold on
for i=21:30
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shift(i)=(index-note_length)/4;
    %smooth
    processed=jc_pinterp(arrayfileD(i).pitches);
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
    plotshiftedD=shifted(i,:)+spacing*(i-20);
    plot(plotshiftedD); xlim([xmin xmax]); ylim([ymin ymax]); title('20-30')
end