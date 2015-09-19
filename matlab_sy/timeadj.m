
%this fxn. moves the start_time and end_time in
%tmvc ahead by the offset amount.
function [tmvc] = timeadj(tmvc, offset, veclist,chopdata)

    %chopdata only shortens the start_time, but does not change the
    %end_time

    %check if time is 7am (.2917)
    %if it is, extend 10 minutes earlier, do not adjust time)
    for ind=1:length(veclist)
        indvl=veclist(ind);
        %start_time
        if(tmvc(indvl,1)-floor(tmvc(indvl,1)))<=.2917
            tmvc(indvl,1)=floor(tmvc(indvl,1))+.25;
        else
            tmvc(indvl,1)=tmvc(indvl,1)+offset;
        end
        
        %end_time
        if(chopdata==0)
            if(tmvc(indvl,2)-floor(tmvc(indvl,2)))>=.875
                tmvc(indvl,2)=floor(tmvc(indvl,2))+.8854
            else
                tmvc(indvl,2)=tmvc(indvl,2)+offset;
            end
        end
    end

%now go through once more and truncate the data if trunc==1
