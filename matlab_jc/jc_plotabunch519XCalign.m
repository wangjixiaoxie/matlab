function jc_plotabunch519XCalign(arrayfileD,avgnote,f)
figure; hold on
spacing=180;

ymin=2200;

ymax=4400; 
figure; hold on
subplot(131); hold on
for i=1:10
    h=xcov(avgnote,f(i).datt);
    [peak,index]=max(h);
    shift=(index-6400)/4;
    starting=300+shift;
    processed=jc_pinterp(arrayfileD(i).pitches(starting:870));
    plotshiftedD=processed+spacing*i;
    plot(plotshiftedD); xlim([0 900]); ylim([ymin ymax]); title('1-10')
end
subplot(132); hold on
for i=10:20
    h=xcov(avgnote,f(i).datt);
    [peak,index]=max(h);
    shift=(index-6400)/4;
    starting=150+shift;
    processed=jc_pinterp(arrayfileD(i).pitches(starting:870));
    plotshiftedD=processed+spacing*(i-10);
    plot(plotshiftedD); xlim([0 900]); ylim([ymin ymax]); title('10-20')
end
subplot(133); hold on
for i=20:30
    h=xcov(avgnote,f(i).datt);
    [peak,index]=max(h);
    shift=(index-6400)/4;
    starting=150+shift;
    processed=jc_pinterp(arrayfileD(i).pitches(starting:870));
    plotshiftedD=processed+spacing*(i-20);
    plot(plotshiftedD); xlim([0 900]); ylim([ymin ymax]); title('20-30')
end