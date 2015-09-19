function shift=jc_plotabunch527(arrayfileD,f)

spacing=180;

ymin=2500;

ymax=4500; 

%Smooth the data
for i=1:length(arrayfileD)
    [holder]=SmoothData(f(i).datt,44100,1);
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
    shift(i)=(index-4096)/4;
    k=1;
    processed=jc_pinterp(arrayfileD(i).pitches);
    for j=1:500
        start=0+shift(i);
        if round(j-start)>0
            proc(i,k)=processed(round(j-start));
        else
            proc(i,k)=3000;
        end
        k=k+1;
    end
    plotshiftedD=proc(i,:)+spacing*i;
    plot(plotshiftedD); xlim([-200 700]); ylim([ymin ymax]); title('1-10')
end

subplot(132); hold on
for i=10:20
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shift(i)=(index-4096)/4;
    starting=180-shift(i);
    processed=jc_pinterp(arrayfileD(i).pitches(starting:700));
    plotshiftedD=processed+spacing*(i-10);
    plot(plotshiftedD); xlim([0 700]); ylim([ymin ymax]); title('10-20')
end
subplot(133); hold on
for i=20:30
    h=xcov(avgnote,smooth(i).smoothed);
    [peak,index]=max(h);
    shift(i)=(index-4096)/4;
    starting=180-shift(i);
    processed=jc_pinterp(arrayfileD(i).pitches(starting:700));
    plotshiftedD=processed+spacing*(i-20);
    plot(plotshiftedD); xlim([0 700]); ylim([ymin ymax]); title('20-30')
end