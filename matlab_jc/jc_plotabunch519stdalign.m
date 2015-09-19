function jc_plotabunch519stdalign(arrayfileD)
figure; hold on
spacing=180;

ymin=2200;

ymax=4400; 
figure; hold on
subplot(131); hold on
for i=1:10
    thresh=jc_getthresh(arrayfileD(i).pitches,100);
    processed=jc_pinterp(arrayfileD(i).pitches(thresh-70:870));
    plotshiftedD=processed+spacing*i;
    plot(plotshiftedD); xlim([0 900]); ylim([ymin ymax]); title('1-10')
end
subplot(132); hold on
for i=10:20
    thresh=jc_getthresh(arrayfileD(i).pitches,100);
    processed=jc_pinterp(arrayfileD(i).pitches(thresh-70:870));
    plotshiftedD=processed+spacing*(i-10);
    plot(plotshiftedD); xlim([0 900]); ylim([ymin ymax]); title('10-20')
end
subplot(133); hold on
for i=20:30
    thresh=jc_getthresh(arrayfileD(i).pitches,100);
    processed=jc_pinterp(arrayfileD(i).pitches(thresh-70:870));
    plotshiftedD=processed+spacing*(i-20);
    plot(plotshiftedD); xlim([0 900]); ylim([ymin ymax]); title('20-30')
end