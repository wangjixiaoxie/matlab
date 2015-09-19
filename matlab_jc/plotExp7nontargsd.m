function plotExp7nontargsd(Experiment7,ravgwin1,window1)
            [timeout,SDout1]=runningstd2(Experiment7.timeACpreCTL(114:end),mean(Experiment7.pitchACpreCTL(window1,114:end)),ravgwin1);
            plot(timeout,SDout1,'Color','k','Linewidth',6)
            [timeout,SDout]=runningstd2(Experiment7.timeAPVwnCTL,mean(Experiment7.pitchAPVwnCTL(window1,:)),ravgwin1);
            plot(timeout,SDout,'Color','r','Linewidth',6)
            [timeout,SDout]=runningstd2(Experiment7.timeAPVCTL,mean(Experiment7.pitchAPVCTL(window1,:)),ravgwin1);
            plot(timeout,SDout,'Color','b','Linewidth',6);
            [timeout,SDout]=runningstd2(Experiment7.timeACpostCTL,mean(Experiment7.pitchACpostCTL(window1,:)),ravgwin1);
            plot(timeout,SDout,'Color','k','Linewidth',6);
            plot([4975 5010],[mean(SDout1) mean(SDout1)],'k')
    
