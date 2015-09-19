

% Example figure
% Positive control learning w/o APV
load /cardinal6/CovertAnalysis/Covert040810.mat
figure;hold on;
i=2
    if length(ExperimentPC(i).time==1)
        window1=ExperimentPC(i).time-60:ExperimentPC(i).time;
    else
        window1=ExperimentPC(i).time;
    end
ravgwin=50;
    plot(ExperimentPC(i).timePre,mean(ExperimentPC(i).pitchPre(window1,:)),'b.')
    plot(ExperimentPC(i).timeWN,mean(ExperimentPC(i).pitchWN(window1,:)),'r.')
     plot(ExperimentPC(i).timePost,mean(ExperimentPC(i).pitchPost(window1,1:length(ExperimentPC(i).timePost))),'k.')    

    plot(runningmedian(ExperimentPC(i).timePre,ravgwin),runningmedian(mean(ExperimentPC(i).pitchPre(window1,:)),ravgwin),'b-','Linewidth',3)
    plot(runningmedian(ExperimentPC(i).timeWN,ravgwin),runningmedian(mean(ExperimentPC(i).pitchWN(window1,:)),ravgwin),'r-','Linewidth',3)
    plot(runningmedian(ExperimentPC(i).timePost,ravgwin),runningmedian(mean(ExperimentPC(i).pitchPost(window1,1:length(ExperimentPC(i).timePost))),ravgwin),'k-','Linewidth',3)    
load /cardinal6/bk76bk63/EPC2.mat
figure;hold on;
window1=
ravgwin=100;
    plot(EPC2.timePre,mean(EPC2.pitchPre(window1,:)),'b.')
    plot(EPC2.timeWN,mean(EPC2.pitchWN(window1,:)),'r.')
     plot(EPC2.timePost,mean(EPC2.pitchPost(window1,1:length(EPC2.timePost))),'k.')    

    plot(runningmedian(EPC2.timePre,ravgwin),runningmedian(mean(EPC2.pitchPre(window1,:)),ravgwin),'g-','Linewidth',3)
    plot(runningmedian(EPC2.timeWN,ravgwin),runningmedian(mean(EPC2.pitchWN(window1,:)),ravgwin),'g-','Linewidth',3)
    plot(runningmedian(EPC2.timePost,20),runningmedian(mean(EPC2.pitchPost(window1,1:length(EPC2.timePost))),20),'g-','Linewidth',3)    



figure;hold on;
for i=1:14
    if length(ExperimentPC(i).time==1)
        window1=ExperimentPC(i).time-30:ExperimentPC(i).time+30;
    else
        window1=ExperimentPC(i).time;
    end
    subplot(4,4,i);hold on;
%     plot(ExperimentPC(i).timePre,mean(ExperimentPC(i).pitchPre(window1,:)),'b.')
%     plot(ExperimentPC(i).timeWN,mean(ExperimentPC(i).pitchWN(window1,:)),'r.')
%     plot(ExperimentPC(i).timePost,mean(ExperimentPC(i).pitchPost(window1,:)),'k.')    

    plot(runningmedian(ExperimentPC(i).timePre,20),runningmedian(mean(ExperimentPC(i).pitchPre(window1,:)),20),'b.')
    plot(runningmedian(ExperimentPC(i).timeWN,20),runningmedian(mean(ExperimentPC(i).pitchWN(window1,:)),20),'r.')
    plot(runningmedian(ExperimentPC(i).timePost,20),runningmedian(mean(ExperimentPC(i).pitchPost(window1,1:length(ExperimentPC(i).timePost))),20),'k.')    
end


ravgwin1=40;       
% Targeted syllable - mean FF
    window1=[140:170];
            % THR ~= 2200Hz
            figure;hold on;subplot(311);hold on;
            plot(Experiment7.timeACpre(114:end),mean(Experiment7.pitchACpre(window1,114:end)),'k.','Markersize',15)
            ratpre=runningmedian(Experiment7.timeACpre(114:end),ravgwin1);
            rafpre=runningmedian(mean(Experiment7.pitchACpre(window1,114:end)),ravgwin1);
            plot(ratpre,rafpre,'Color','g','Linewidth',6)
            plot(Experiment7.timeAPV,mean(Experiment7.pitchAPV(window1,:)),'b.','Markersize',15)
            plot(Experiment7.timeAPVwn,mean(Experiment7.pitchAPVwn(window1,:)),'r.','Markersize',15)
            ratapvwn=runningmedian(Experiment7.timeAPVwn,ravgwin1);
            rafapvwn=runningmedian(mean(Experiment7.pitchAPVwn(window1,:)),ravgwin1);
            plot(ratapvwn,rafapvwn,'g','Linewidth',6)
            ratapv=runningmedian(Experiment7.timeAPV,ravgwin1);
            rafapv=runningmedian(mean(Experiment7.pitchAPV(window1,:)),ravgwin1);
            plot(ratapv,rafapv,'g','Linewidth',6)
            plot(Experiment7.timeACpost,(mean(Experiment7.pitchACpost(window1,:))),'k.','Markersize',15)
            plot(Experiment7.timeACpost(1:end),(mean(Experiment7.pitchACpost(window1,1:end))),'k.','Markersize',15)
            ratpost=runningmedian(Experiment7.timeACpost(1:end),ravgwin1);
            rafpost=runningmedian(mean(Experiment7.pitchACpost(window1,1:end)),ravgwin1);
            ratpost=runningmedian(Experiment7.timeACpost,ravgwin1);
            rafpost=runningmedian(mean(Experiment7.pitchACpost(window1,:)),ravgwin1);
            plot(ratpost,rafpost,'g','Linewidth',6)
            %avgescapes=mean(mean(pitchE(window1,:)))-mean(mean(Experiment7.pitchACpre(window1,:)));
            plot([4975 5010],[mean(rafpre) mean(rafpre)],'k')
            plot([min(ratapvwn) max(ratapvwn)],[2200 2200],'k')
             ylim([2000 2400])           
 % Targeted syllable - s.d. of FF           
    subplot(312);hold on;
            ratpre=runningmedian(Experiment7.timeACpre(114:end),ravgwin1);
            rafpre=runningiqr(mean(Experiment7.pitchACpre(window1,114:end)),ravgwin1);
            plot(ratpre,rafpre,'Color','k','Linewidth',6)
            ratapvwn=runningmedian(Experiment7.timeAPVwn,ravgwin1);
            rafapvwn=runningiqr(mean(Experiment7.pitchAPVwn(window1,:)),ravgwin1);
            plot(ratapvwn,rafapvwn,'r','Linewidth',6)
            ratapv=runningmedian(Experiment7.timeAPV,ravgwin1);
            rafapv=runningiqr(mean(Experiment7.pitchAPV(window1,:)),ravgwin1);
            plot(ratapv,rafapv,'b','Linewidth',6)
            ratpost=runningmedian(Experiment7.timeACpost(1:end),ravgwin1);
            rafpost=runningiqr(mean(Experiment7.pitchACpost(window1,1:end)),ravgwin1);
            ratpost=runningmedian(Experiment7.timeACpost,ravgwin1);
            rafpost=runningiqr(mean(Experiment7.pitchACpost(window1,:)),ravgwin1);
            plot(ratpost,rafpost,'k','Linewidth',6)
            %avgescapes=mean(mean(pitchE(window1,:)))-mean(mean(Experiment7.pitchACpre(window1,:)));
            plot([4975 5010],[mean(rafpre) mean(rafpre)],'k')
            %plot([min(ratapvwn) max(ratapvwn)],[2200 2200],'k')           

 % Non-targeted syllable - mean FF           
            window1=[180:340];
            subplot(313);hold on;
            plot(Experiment7.timeACpreCTL(114:end),mean(Experiment7.pitchACpreCTL(window1,114:end)),'k.','Markersize',15)
            ratpre=runningmedian(Experiment7.timeACpreCTL(114:end),ravgwin1);
            rafpre=runningmedian(mean(Experiment7.pitchACpreCTL(window1,114:end)),ravgwin1);
            plot(ratpre,rafpre,'Color','g','Linewidth',6)
            plot(Experiment7.timeAPVCTL,mean(Experiment7.pitchAPVCTL(window1,:)),'b.','Markersize',15)
            plot(Experiment7.timeAPVwnCTL,mean(Experiment7.pitchAPVwnCTL(window1,:)),'r.','Markersize',15)
            ratapvwn=runningmedian(Experiment7.timeAPVwnCTL,ravgwin1);
            rafapvwn=runningmedian(mean(Experiment7.pitchAPVwnCTL(window1,:)),ravgwin1);
            plot(ratapvwn,rafapvwn,'g','Linewidth',6)
            ratapv=runningmedian(Experiment7.timeAPVCTL,ravgwin1);
            rafapv=runningmedian(mean(Experiment7.pitchAPVCTL(window1,:)),ravgwin1);
            plot(ratapv,rafapv,'g','Linewidth',6)
            plot(Experiment7.timeACpostCTL(1:end),(mean(Experiment7.pitchACpostCTL(window1,1:end))),'k.','Markersize',15)
            %plot(Experiment7.timeACpostCTL,(mean(Experiment7.pitchACpostCTL(window1,:))),'k.','Markersize',15)
            ratpost=runningmedian(Experiment7.timeACpostCTL(1:end),ravgwin1);
            rafpost=runningmedian(mean(Experiment7.pitchACpostCTL(window1,1:end)),ravgwin1);
            plot(ratpost,rafpost,'g','Linewidth',6)
            %avgescapes=mean(mean(pitchE(window1,:)))-mean(mean(Experiment7.pitchACpre(window1,:)));
            plot([4975 5010],[mean(mean(Experiment7.pitchACpreCTL(window1,114:end))) mean(mean(Experiment7.pitchACpreCTL(window1,114:end)))],'k')
            ylim([1900 2300])

% Non-targeted syllable - s.d. of FF
%           subplot(313);hold  on
%             ravgwin1=80;
%             window1=[180:340];
% %             plot(Experiment7.timeACpreCTL(114:end),mean(Experiment7.pitchACpreCTL(window1,114:end)),'k.','Markersize',15)
%             ratpre=runningmedian(Experiment7.timeACpreCTL(114:end),ravgwin1);
%             rafpre=runningiqr(mean(Experiment7.pitchACpreCTL(window1,114:end)),ravgwin1);
%             plot(ratpre,rafpre,'Color','k','Linewidth',6)
% %             plot(Experiment7.timeAPVCTL,mean(Experiment7.pitchAPVCTL(window1,:)),'b.','Markersize',15)
% %             plot(Experiment7.timeAPVwnCTL,mean(Experiment7.pitchAPVwnCTL(window1,:)),'r.','Markersize',15)
%             ratapvwn=runningmedian(Experiment7.timeAPVwnCTL,ravgwin1);
%             rafapvwn=runningiqr(mean(Experiment7.pitchAPVwnCTL(window1,:)),ravgwin1);
%             plot(ratapvwn,rafapvwn,'r','Linewidth',6)
%             ratapv=runningmedian(Experiment7.timeAPVCTL,ravgwin1);
%             rafapv=runningiqr(mean(Experiment7.pitchAPVCTL(window1,:)),ravgwin1);
%             plot(ratapv,rafapv,'g','Linewidth',6)
% %             plot(Experiment7.timeACpostCTL(1:end),(mean(Experiment7.pitchACpostCTL(window1,1:end))),'k.','Markersize',15)
%             %plot(Experiment7.timeACpostCTL,(mean(Experiment7.pitchACpostCTL(window1,:))),'k.','Markersize',15)
%             ratpost=runningmedian(Experiment7.timeACpostCTL(1:end),ravgwin1);
%             rafpost=runningiqr(mean(Experiment7.pitchACpostCTL(window1,1:end)),ravgwin1);
%             plot(ratpost,rafpost,'b','Linewidth',6)
%             %avgescapes=mean(mean(pitchE(window1,:)))-mean(mean(Experiment7.pitchACpre(window1,:)));
%             %plot([4975 5010],[mean(mean(Experiment7.pitchACpreCTL(window1,114:end))) mean(mean(Experiment7.pitchACpreCTL(window1,114:end)))],'k')
