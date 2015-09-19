function plotExp7targsd(Experiment7,ravgwin1,window1)
            [timeout1,SDout1]=runningstd2(Experiment7.timeACpre(114:end),mean(Experiment7.pitchACpre(window1,114:end)),ravgwin1);
            plot(timeout1,SDout1,'Color','k','Linewidth',6)
            [timeout,SDout]=runningstd2(Experiment7.timeAPVwn,mean(Experiment7.pitchAPVwn(window1,:)),ravgwin1);
            plot(timeout,SDout,'Color','r','Linewidth',6)
            [timeout,SDout]=runningstd2(Experiment7.timeAPV,mean(Experiment7.pitchAPV(window1,:)),ravgwin1);
            plot(timeout,SDout,'Color','b','Linewidth',6);
            [timeout,SDout]=runningstd2(Experiment7.timeACpost,mean(Experiment7.pitchACpost(window1,:)),ravgwin1);
            plot(timeout,SDout,'Color','k','Linewidth',6);
            plot([4975 5010],[mean(SDout1) mean(SDout1)],'k')
            %plot([min(ratapvwn) max(ratapvwn)],[2200 2200],'k')           
