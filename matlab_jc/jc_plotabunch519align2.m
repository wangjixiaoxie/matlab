function jc_plotabunch519align2(arrayfileD)

spacing=180;

ymin=2200;

ymax=4400; 

for i=1:length(arrayfileD(1).pitches)
    sumnote(i)=0;
    for j=1:length(arrayfileD)
        sumnote(i)=sumnote(i)+arrayfileD(j).pitches(i);
    end
    avgnote(i)=sumnote(i)/length(arrayfileD);
end
avgnote=arrayfileD(1).pitches;
figure; hold on
subplot(131); hold on
for i=1:10
    h=xcov(avgnote,arrayfileD(i).pitches);
    [peak,index]=max(h);
    shift=(index-870);
    starting=250-shift;
    processed=jc_pinterp(arrayfileD(i).pitches(starting:870));
    plotshiftedD=processed+spacing*i;
    plot(plotshiftedD); xlim([0 900]); ylim([ymin ymax]); title('1-10')
end
subplot(132); hold on
for i=10:20
    h=xcov(avgnote,arrayfileD(i).pitches);
    [peak,index]=max(h);
    shift=(index-870);
    starting=100+shift;
    processed=jc_pinterp(arrayfileD(i).pitches(starting:870));
    plotshiftedD=processed+spacing*(i-10);
    plot(plotshiftedD); xlim([0 900]); ylim([ymin ymax]); title('10-20')
end
subplot(133); hold on
for i=20:30
    h=xcov(avgnote,arrayfileD(i).pitches);
    [peak,index]=max(h);
    shift=(index-870);
    starting=100+shift;
    processed=jc_pinterp(arrayfileD(i).pitches(starting:870));
    plotshiftedD=processed+spacing*(i-20);
    plot(plotshiftedD); xlim([0 900]); ylim([ymin ymax]); title('20-30')
end