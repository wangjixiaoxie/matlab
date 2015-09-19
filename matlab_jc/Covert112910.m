
% Experiment 12
load /cardinal6/CovertAnalysis/Covert040810.mat
load /cardinal6/CovertAnalysis/Example12_112910.mat
figure;hold on;
subplot(121);hold on;
% equal sampling
indACpre=1:5:length(fvalsACpre);
indAPV=1:5:length(fvalsAPV);
indAPV=indAPV(find(indAPV>256));
indAPVwn=1:5:length(fvalsAPVwn);
indACpost=1:5:length(fvalsACpost);
indACpost=indACpost(find(indACpost>113));
%
mnbase=mean(mean(pitchACpre(400:350,indACpre)));
plot(timing4(fvalsACpre(indACpre)),mean(pitchACpre(400:350,indACpre))-mnbase,'k.','Markersize',15)
plot(timing4(fvalsAPV(indAPV)),mean(pitchAPV(400:350,indAPV))-mnbase,'b.','Markersize',15) % exclude first hour
plot(timing4(fvalsAPVwn(indAPVwn)),mean(pitchAPVwn(400:350,indAPVwn))-mnbase,'r.','Markersize',15)
plot(timing4(fvalsACpost(indACpost)),mean(pitchACpost(400:350,indACpost))-mnbase,'k.','Markersize',15) % next morning
plot([min(timing4(fvalsAPVwn)) max(timing4(fvalsAPVwn))],[2610-mnbase 2610-mnbase],'b')
plot([min(timing4(fvalsACpre))-5 max(timing4(fvalsACpost))+5],[0 0],'k')
mnpost=mean(mean(pitchACpost(400:350,114:end)))-mnbase;
plot([6640 6636],[mnpost mnpost],'k')
ylim([-150 150])
xlim([6605 6636])


% % Experiment 13
% load /cardinal6/CovertAnalysis/Example12.mat
% figure;hold on;
% plot(Example12.timeACpre,mean(Example12.pitchACpre(400:350,:)),'k.','Markersize',15)
% plot(Example12.timeAPV,mean(Example12.pitchAPV(400:350,:)),'b.','Markersize',15) % exclude first hour
% plot(Example12.timeAPVwn,mean(Example12.pitchAPVwn(400:350,:)),'r.','Markersize',15)
% plot(Example12.timeACpost,mean(Example12.pitchACpost(400:350,:)),'k.','Markersize',15) % next morning
% plot([min(Example12.timeAPVwn) max(Example12.timeAPVwn)],[2610 2610],'k')
% mnbase=mean(mean(Example12.pitchACpre(400:350,:)));
% plot([min(Example12.timeACpre) 6525],[mnbase mnbase],'k')
% mnpost=mean(mean(Example12.pitchACpost(400:350,36:end)));
% plot([6510 6525],[mnpost mnpost],'k')
% ylim([2440 2740])
% xlim([6594 6636])

load /cardinal6/CovertAnalysis/Experiment13.mat
    window11=[150:340];
    timePre=timing4(Experiment13.fvACpre);
    timePost=timing4(Experiment13.fvACpost(114:end));
    mpre=mean(mean(Experiment13.pitchACpre(window11,:)));
    mpost=mean(mean(Experiment13.pitchACpost(window11,114:end)));
figure;hold on;
subplot(221);hold on;
    plot(timing4(Experiment13.fvACpre),mean(Experiment13.pitchACpre(window11,:)),'k.','Markersize',15)
    plot(runningaverage(timing4(Experiment13.fvACpre),50),runningaverage(mean(Experiment13.pitchACpre(window11,:)),50),'Color','r','Linewidth',6) % next morning
    plot(timing4(Experiment13.fvAPV(143:end)),mean(Experiment13.pitchAPV(window11,143:end)),'b.','Markersize',15) % exclude first hour
    plot(timing4(Experiment13.fvAPVwn),mean(Experiment13.pitchAPVwn(window11,:)),'r.','Markersize',15)
    plot(timing4(Experiment13.fvACpost(114:end)),mean(Experiment13.pitchACpost(window11,114:end)),'k.','Markersize',15) % next morning
    plot(runningaverage(timing4(Experiment13.fvACpost((114:end))),50),runningaverage(mean(Experiment13.pitchACpost(window11,(114:end))),50),'Color','r','Linewidth',6) % next morning
    plot([min(timePre) max(timePost)],[mpre mpre],'k')
    plot([min(timePost) max(timePost)],[mpost mpost],'r','Linewidth',3)
    ylim([2350 2750])
subplot(223);hold on;
    window12=[1100 1400];
    timePre=timing4(Experiment13.fvACpre);
    timePost=timing4(Experiment13.fvACpost(114:end));
    mpre=median(mean(Experiment13.pitchACpre(window12,:)));
    mpost=median(mean(Experiment13.pitchACpost(window12,114:end)));
    indACpre=find(std(Experiment13.pitchACpre(window12,:))<100);
    indAPV=find(std(Experiment13.pitchAPV(window12,:))<100);        
    indACpost=find(std(Experiment13.pitchACpost(window12,:))<100);    
    plot(timing4(Experiment13.fvACpre(indACpre)),mean(Experiment13.pitchACpre(window12,indACpre)),'k.','Markersize',15)
    plot(runningaverage(timing4(Experiment13.fvACpre(indACpre)),50),runningaverage(mean(Experiment13.pitchACpre(window12,indACpre)),50),'Color','r','Linewidth',6) % next morning
    plot(timing4(Experiment13.fvAPV(indAPV(143:end))),mean(Experiment13.pitchAPV(window12,indAPV(143:end))),'b.','Markersize',15) % exclude first hour
    plot(timing4(Experiment13.fvAPVwn),mean(Experiment13.pitchAPVwn(window12,:)),'r.','Markersize',15)
    plot(timing4(Experiment13.fvACpost(indACpost(113:end))),mean(Experiment13.pitchACpost(window12,indACpost(113:end))),'k.','Markersize',15) % next morning
    plot(runningaverage(timing4(Experiment13.fvACpost(indACpost(113:end))),50),runningaverage(mean(Experiment13.pitchACpost(window12,indACpost(113:end))),50),'Color','r','Linewidth',6) % next morning
    plot([min(timePre) max(timePost)],[mpre mpre],'k')
    plot([min(timePost) max(timePost)],[mpost mpost],'r','Linewidth',3)
    ylim([2350 2750])
    
% Experiment 14
window1=[400:340];
figure;subplot(222);hold on;
plot(Experiment(14).timeACpre,mean(Experiment(14).pitchACpre(window1,:)),'k.','Markersize',15)
plot(runningaverage(Experiment(14).timeACpre,50),runningaverage(mean(Experiment(14).pitchACpre(window1,:)),50),'Color','r','Linewidth',6)
plot(Experiment(14).timeAPV(11:end),mean(Experiment(14).pitchAPV(window1,11:end)),'b.','Markersize',15)
plot(Experiment(14).timeAPVwn,mean(Experiment(14).pitchAPVwn(window1,:)),'r.','Markersize',15)
plot(Experiment(14).timeACpost(250:382),(mean(Experiment(14).pitchACpost(window1,250:382))),'k.','Markersize',15)
plot(runningaverage(Experiment(14).timeACpost(250:382),100),runningaverage(mean(Experiment(14).pitchACpost(window1,250:382)),100),'Color','r','Linewidth',6)
%avgescapes=mean(mean(pitchE(window1,:)))-mean(mean(Experiment(14).pitchACpre(window1,:)));
plot([8165 8180],[mean(mean(Experiment(14).pitchACpost(window1,250:382))) mean(mean(Experiment(14).pitchACpost(window1,250:382)))] ,'Color','r','Linewidth',3)
plot([8148 8155],[avgescapes avgescapes])
xlim([8140 8180])
ylim([2950 3400])    
plot([8140 8180],[mean(mean(Experiment(14).pitchACpre(window1,:))) mean(mean(Experiment(14).pitchACpre(window1,:)))],'k')
subplot(224);hold on;
plot(Experiment(14).timeACpreCTL,(mean(Experiment(14).pitchACpreCTL(window1,:))),'k.','Markersize',15)
%plot(Experiment(14).timeAPVCTL(11:end),(mean(Experiment(14).pitchAPVCTL(window1,11:end))-mean(mean(Experiment(14).pitchACpreCTL(window1,:)))),'b.','Markersize',15)
%plot(Experiment(14).timeAPVwnCTL,(mean(Experiment(14).pitchAPVwnCTL(window1,:))-mean(mean(Experiment(14).pitchACpreCTL(window1,:)))),'r.','Markersize',15)
plot(Experiment(14).timeACpostCTL(250:382),(mean(Experiment(14).pitchACpostCTL(window1,250:382))),'k.','Markersize',15)
plot(runningaverage(Experiment(14).timeACpostCTL(250:382),50),runningaverage(mean(Experiment(14).pitchACpostCTL(window1,250:382)),50),'Color','r','Linewidth',6)
%avgescapes=mean(mean(pitchE(window1,:)))-mean(mean(Experiment(14).pitchACpre(window1,:)));
plot([8165 8180],[mean(mean(Experiment(14).pitchACpostCTL(window1,250:382))) mean(mean(Experiment(14).pitchACpostCTL(window1,250:382)))] ,'Color','r','Linewidth',3)
%plot([8148 8155],[avgescapes avgescapes])
xlim([8140 8180])
plot([8140 8180],[mean(mean(Experiment(14).pitchACpreCTL(window1,:))) mean(mean(Experiment(14).pitchACpreCTL(window1,:)))],'Color','k')    
ylim([2950 3400])    

allpts=[mean(Experiment7.pitchACpre(window1,114:end)) mean(Experiment7.pitchAPV(window1,:)) ...
    mean(Experiment7.pitchAPVwn(window1,:)) mean(Experiment7.pitchACpost(window1,1:end))]; 
allpts2=[zeros(1,50) allpts-mean(mean(Experiment7.pitchACpre(window1,114:end)))];
