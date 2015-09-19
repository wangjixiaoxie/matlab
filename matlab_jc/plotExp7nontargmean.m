function plotExp7nontargmean(Experiment7,ravgwin1,window1)
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