% False positives part 1 - mean autocorr from white noise
[acS,widthsS,xpsS,psdS]=jc_flucWN;
[acI,widthsI,xpsI,psdI]=jc_flucWNifd;
figure;plot(acI')
hold on;plot(acS(1,:),'k')
hold on;plot(acS(2,:),'k')
hold on;plot(acS(3,:),'k')

% False positives part 2 - distn of autocorr widths from white noise
nbins=0:1:30;
figure;
subplot(611)
hist(widthsI(1,:)/8,nbins)
xlim([0 30])
ylim([0 30])
subplot(612)
hist(widthsI(2,:)/8,nbins)
xlim([0 30])
ylim([0 30])
subplot(613)
hist(widthsI(3,:)/8,nbins)
xlim([0 30])
ylim([0 30])

subplot(614)
hist(widthsS(1,:)/8,nbins)
xlim([0 30])
ylim([0 30])
subplot(615)
hist(widthsS(2,:)/8,nbins)
xlim([0 30])
ylim([0 30])
subplot(616)
hist(widthsS(3,:)/8,nbins)
xlim([0 30])
ylim([0 30])

% Part 3 - PSD contributions
figure;plot(xpsS,psdS);hold on;plot(xpsI,psdI,'k')

% Part 4 - Consent Plots
for k=1:11
    cpZF(k).data=ConsentPlot(AlldataZFlesion(k).rawdataUDpre(1,:),32000,1024,1020,1,3,5,5);
end

FS=[32000 44100 44100 44100 44100 32000 32000 32000 32000];
for j=1:9
    cpBF(j).data=ConsentPlot(AlldataBFlesion(j).rawdataUDpre(1,:),FS(j),1024,1020,1,3,5,5);
end