%screen script
%to look for candidate birdind

%not including r34w20 yet- add later
ct=1;
for bsind=1:4
    crbs=sumbs(bsind);
    for runind=1:length(crbs.ptvls);
        crstmn=mean(crbs.STIMTIME{runind});
        crststd=std(crbs.STIMTIME{runind});
        crdel=crbs.del{runind};
        
        if( (crstmn-2*crststd+crdel+.08)<.11)
            screenout(ct).bsind=bsind;
            screenout(ct).runvl=runind;
            ct=ct+1;
        end
    end
end
  