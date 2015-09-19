
%this fxn. moves the start_time and end_time in
%tmvc ahead by the offset amount.
function [tmvc] = timeadj(tmvc, offset, veclist,chopdata)
    aclist=avls.aclist
    mulist=avls.mulist;
    acoff=avls.acoffset;
    muoff=avls.muoffset;
    tmvc=avls.rawtimes;
    %check if time is 7am (.2917)
    %if it is, extend 10 minutes earlier, do not adjust time)
    for ind=1:length(mulist)
        muindvl=mulist(ind);
        preindvl=aclist(ind,1);
        pstindvl=aclist(ind,2);
        %start_time
        if(tmvc(preindvl,1)-floor(tmvc(preindvl,1)))<=.2917
            tmvcout(preindvl,1)=floor(tmvc(preindvl,1)+.2812);
        elseif(tmvc(pstindvl,2)-floor(tmvc(pstindvl,2)))>=.875;
            tmvcout(pstindvl,2)=floor(tmvc(pstindvl,2))+.8854;
            end
        
        tmvcout(muindvl,1)=tmvc(muindvl,1)+muoff;
        tmvcout(pstindvl,1)=tmvc(pstindvl,1)+acoff;
    end
        

%now go through once more and truncate the data if trunc==1
