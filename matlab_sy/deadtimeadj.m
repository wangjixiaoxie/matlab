
%this fxn. moves the start_time and end_time in
%tmvc ahead by the offset amount.
function [tmvc] = deadtimeadj(avls)
    aclist=avls.aclist
    mulist=avls.mulist
    deadtm=avls.deadtm
    tmvc=avls.rawtimes
    %chopdata only shortens the start_time, but does not change the
    %end_time

    %check if time is 7am (.2917)
    %if it is, extend 10 minutes earlier, do not adjust time)
    muruns=find(mulist>0);
    for ind=1:length(muruns)
        indvl=muruns(ind);
        muind=mulist(indvl);
        acpreind=aclist(indvl,1);
        acpstind=aclist(indvl,2);
        %mv acpretime back by deadtime;
        tmvc(acpreind,2)=tmvc(acpreind,2)+deadtm;
        %mv acpst back by deadtm
        tmvc(acpstind,1)=tmvc(acpstind,1)+deadtm;
        %mv mupre back
        tmvc(muind,2)=tmvc(muind,2)+deadtm;
        tmvc(muind,1)=tmvc(muind,1)+deadtm;
    end
        
        