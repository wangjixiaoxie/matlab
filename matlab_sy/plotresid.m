
    clear ps
    outstrinds=[3 4 5 6  ]
    subplotinds=[4 7 10 13 ]
    stimdiff=[.045 .035 .025 .015]
    for ii=1:length(outstrinds)
        ctcontours=outstr.combctcon{outstrinds(ii)};
        fbcontours=outstr.combfbcon{outstrinds(ii)};
        axcon(ii)=subplot(5,3,subplotinds(ii))
        ps.ax=gca();
        ps.basind=1
        ps.TIMEBNDS=[.07 .11]
        
        ps.plotstim=1;
%         ps.stimvals=datnew;
        ps.stimdiff=stimdiff(ii);
        ps.STIMGAIN=.000001;
        ps.STIMOFFSET=2.3
        plotcontour2(ctcontours,fbcontours,avls,ps)
        axis square
    end
    linkaxes(axcon);
    axis([0.02 0.12 2.25 2.4])
    
