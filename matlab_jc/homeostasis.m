% bk72bk64 --- look after 1 month - august 8th - postlesion2/trigtest2 -
% huge circadian changes

figure;hold on;
plot(runningaverage(timing3(fvalsApre1(ind1)),20),runningaverage(mean(pitchApre1(1000:1300,ind1)),20),'*')
plot(runningaverage(timing3(fvalsApre2),20),runningaverage(mean(pitchApre2(1000:1300,:)),20),'*')
plot(runningaverage(timing3(fvalsApost1),20),runningaverage(mean(pitchApost1(1000:1300,:)),20),'*')
plot(runningaverage(timing3(fvalsApost2),20),runningaverage(mean(pitchApost2(1000:1300,:)),20),'*')
plot(runningaverage(timing3(fvalsApostB1),5),runningaverage(mean(pitchApostB1(1000:1300,:)),5),'*')
plot(runningaverage(timing3(fvalsApostB2),5),runningaverage(mean(pitchApostB2(1000:1300,:)),5),'*')

