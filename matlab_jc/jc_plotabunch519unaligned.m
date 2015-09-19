function jc_plotabunch519unaligned(arrayfileD)
figure; hold on
spacing=180;

ymin=3000;

ymax=6000; 
figure; hold on
subplot(131); hold on
for i=1:10
    %thresh=jc_getthresh(arrayfileD(i).pitches,100);
    processed=jc_pinterp(arrayfileD(i).pitches(1:645));
    plotshiftedD=processed+spacing*i;
    plot(plotshiftedD); xlim([0 700]); ylim([ymin ymax]); title('1-10')
end
subplot(132); hold on
for i=10:20
   % thresh=jc_getthresh(arrayfileD(i).pitches,100);
    processed=jc_pinterp(arrayfileD(i).pitches(1:645));
    plotshiftedD=processed+spacing*(i-10);
    plot(plotshiftedD); xlim([0 700]); ylim([ymin ymax]); title('10-20')
end
subplot(133); hold on
for i=20:30
   % thresh=jc_getthresh(arrayfileD(i).pitches,100);
    processed=jc_pinterp(arrayfileD(i).pitches(1:645));
    plotshiftedD=processed+spacing*(i-20);
    plot(plotshiftedD); xlim([0 700]); ylim([ymin ymax]); title('20-30')
end