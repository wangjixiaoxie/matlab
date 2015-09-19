function jc_plotabunch512(arrayfileD)
figure; hold on
spacing=180;
%sample_min=1;
%sample_max=25;
ymin=2200;
%bulge=400; % Fits stuff on graph
ymax=4400; %ymin+bulge+spacing*(sample_max-sample_min);
figure; hold on
subplot(131); hold on
for i=1:10
    thresh=jc_getthresh(arrayfileD(i).pitches,100);
    processed=jc_pinterp(arrayfileD(i).pitches(thresh:640));
    plotshiftedD=processed+spacing*i;
    plot(plotshiftedD); xlim([0 800]); ylim([ymin ymax]); title('1-10')
end
subplot(132); hold on
for i=10:20
    thresh=jc_getthresh(arrayfileD(i).pitches,100);
    processed=jc_pinterp(arrayfileD(i).pitches(thresh:640));
    plotshiftedD=processed+spacing*(i-10);
    plot(plotshiftedD); xlim([0 800]); ylim([ymin ymax]); title('10-20')
end
subplot(133); hold on
for i=20:30
    thresh=jc_getthresh(arrayfileD(i).pitches,100);
    processed=jc_pinterp(arrayfileD(i).pitches(thresh:640));
    plotshiftedD=processed+spacing*(i-20);
    plot(plotshiftedD); xlim([0 800]); ylim([ymin ymax]); title('20-30')
end