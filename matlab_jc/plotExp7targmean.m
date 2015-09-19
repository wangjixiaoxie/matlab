function plotExp7targmean(Experiment7,ravgwin1,window1)
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