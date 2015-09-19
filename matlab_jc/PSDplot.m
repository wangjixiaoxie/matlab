% load /cardinal4/PSDdata.mat
% normBiFL(1).data --- ACSF
% normBiFL(1).dataC --- INA
% normBFL(1).data --- UDpre
% normBFL(1).dataB --- Dpre
% normBFL(1).dataC --- UDpost
% Each of the files uses 260pts centered at the middle of the note because
% the objective is to prevent contributions of crappy stuff at the edges.
% The sole exception is #4 from the normBiFL inactivation data set because
% this note is extremely short.

%%% This figure shows the raw data
% powerBFvsZFrawdata.fig
figure;hold on;
for i=1:4
    hold on;plot(x,psdBFL(i).data./psdBFL(i).dataC,'r')
end
for i=1:3
    hold on;plot(x,psdBiFL(i).data./psdBiFL(i).dataC,'g')
end
hold on;plot(x(1:256),psdBiFL(4).data./psdBiFL(4).dataC,'g')
for i=1:6
    hold on;plot(x,psdZFL(i).data./psdZFL(i).dataC,'b')
end

%%% Power at 32ms divided by power at 16ms
    % powerBFvsZFsummary.fig
figure;hold on
pt1=find(x==32);
pt2=find(x==16);
for i=1:4
    BiRatio(i).data=psdBiFL(i).data./psdBiFL(i).dataC;
    BLRatio(i).data=psdBFL(i).data./psdBFL(i).dataC;
end
for i=1:6
    ZLRatio(i).data=psdZFL(i).data./psdZFL(i).dataC;
end
for i=1:4
    plot(1,BiRatio(i).data(pt1)./BiRatio(i).data(pt2),'+')
    plot(2,BLRatio(i).data(pt1)./BLRatio(i).data(pt2),'+')
end
for i=1:6
    plot(3,ZLRatio(i).data(pt1)./ZLRatio(i).data(pt2),'+')
end
    % plot a line at unity
    for i=1:1000
        zz(i)=i/200;
    end
    hold on;plot(zz,1)



% normBFL(1).data=residuals(AlldataBFlesion(1).pitchUDpre(570:830,1:50));
% normBFL(1).dataB=residuals(AlldataBFlesion(1).pitchDpre(570:830,1:50));
% normBFL(1).dataC=residuals(AlldataBFlesion(1).pitchUDpost(570:830,1:50));
% normBFL(2).data=residuals(AlldataBFlesion(1).pitchUDpre(320:580,1:45));
% normBFL(2).data=residuals(AlldataBFlesion(2).pitchUDpre(320:580,1:45));
% normBFL(2).dataB=residuals(AlldataBFlesion(2).pitchDpre(320:580,1:45));
% normBFL(2).dataC=residuals(AlldataBFlesion(2).pitchUDpost(320:580,1:45));
% normBFL(3).data=residuals(AlldataBFlesion(3).pitchUDpre(420:680,1:30));
% normBFL(3).dataB=residuals(AlldataBFlesion(3).pitchDpre(420:680,1:30));
% normBFL(3).dataC=residuals(AlldataBFlesion(3).pitchUDpost(420:680,1:30));
% normBFL(4).data=residuals(AlldataBFlesion(4).pitchUDpre(300:560,1:22));
% normBFL(4).dataB=residuals(AlldataBFlesion(4).pitchDpre(300:560,1:22));
% normBFL(4).dataC=residuals(AlldataBFlesion(4).pitchUDpost(300:560,1:22));
% normZFL(1).data=residuals(AlldataZFlesion(1).pitchUDpre(150:410,1:50));
% normZFL(1).dataB=residuals(AlldataZFlesion(1).pitchDpre(150:410,1:50));
% normZFL(1).dataC=residuals(AlldataZFlesion(1).pitchUDpost(150:410,1:50));
% normZFL(2).data=residuals(AlldataZFlesion(1).pitchUDpre(700:960,1:45));
% normZFL(2).dataB=residuals(AlldataZFlesion(1).pitchDpre(700:960,1:45));
% normZFL(2).dataC=residuals(AlldataZFlesion(1).pitchUDpost(700:960,1:45));
% normZFL(2).data=residuals(AlldataZFlesion(2).pitchUDpre(700:960,1:45));
% normZFL(2).dataB=residuals(AlldataZFlesion(2).pitchDpre(700:960,1:45));
% normZFL(2).dataC=residuals(AlldataZFlesion(2).pitchUDpost(700:960,1:45));
% normZFL(3).data=residuals(AlldataZFlesion(3).pitchUDpre(510:770,1:50));
% normZFL(3).dataB=residuals(AlldataZFlesion(3).pitchDpre(510:770,1:50));
% normZFL(3).dataC=residuals(AlldataZFlesion(3).pitchUDpost(510:770,1:50));
% normZFL(4).data=residuals(AlldataZFlesion(3).pitchUDpre(680:940,1:44));
% normZFL(4).dataB=residuals(AlldataZFlesion(4).pitchDpre(680:940,1:44));
% normZFL(4).data=residuals(AlldataZFlesion(4).pitchUDpre(680:940,1:44));
% normZFL(4).dataC=residuals(AlldataZFlesion(4).pitchUDpost(680:940,1:44));
% normZFL(5).data=residuals(AlldataZFlesion(5).pitchUDpre(150:410,1:55));
% normZFL(5).dataB=residuals(AlldataZFlesion(5).pitchDpre(150:410,1:55));
% normZFL(5).dataC=residuals(AlldataZFlesion(5).pitchUDpost(150:410,1:55));
% normZFL(6).data=residuals(AlldataZFlesion(6).pitchUDpre(450:710,1:55));
% normZFL(6).dataB=residuals(AlldataZFlesion(6).pitchDpre(450:710,1:55));
% normZFL(6).dataC=residuals(AlldataZFlesion(6).pitchUDpost(450:710,1:55));
% 
% normBiFL(1).dataC=residuals(Goodinactivations.INA(2).pitchSTFT(500:760,1:12));
% normBiFL(1).data=residuals(Goodinactivations.ACSF(2).pitchSTFT(500:760,1:12));
% normBiFL(2).data=residuals(Goodinactivations.ACSF(3).pitchSTFT(300:560,1:33));
% normBiFL(2).dataC=residuals(Goodinactivations.INA(3).pitchSTFT(300:560,1:33));
% normBiFL(3).data=residuals(Goodinactivations.ACSF(4).pitchSTFT(750:1010,1:100));
% normBiFL(3).dataC=residuals(Goodinactivations.INA(4).pitchSTFT(750:1010,1:100));
% normBiFL(4).data=residuals(Goodinactivations.ACSF(4).pitchSTFT(280:410,1:60));
% normBiFL(4).dataC=residuals(Goodinactivations.INA(4).pitchSTFT(280:410,1:60));
% 
% [x,psdBiFL(i).data]=jc_psd(normBiFL(i).data);
% 
% for i=1:4
% hold on;plot(x,psdBiFL(i).data./psdBiFL(i).dataC,'g')
% end