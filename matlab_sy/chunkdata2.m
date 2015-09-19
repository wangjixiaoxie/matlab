function [avls]=chunkdata2(avls,graphvals,useofst);
    %start by taking avls.vls{ntind} and finding the indices that match
    %these time periods
    %calls fxn sorttimes, getmean, getstd, and getvals.
    for ntind=1:length(avls.NT)
        
        avls.vls{ntind}=getvals(avls.fvcomb{ntind},1,'TRIG')
        
        avls.rawvls{ntind}=sort_times(avls.vls{ntind},avls.rawtimes);
        avls.allvls{ntind}=sort_times(avls.vls{ntind},avls.alltimes)
        avls.adjvls{ntind}=sort_times(avls.vls{ntind},avls.adjtimes);
%         if(~isempty(avls.diron))
%             avls.dirvls{ntind}=getvals(avls.fvcombdir{ntind},1,'TRIG');
%           
%             avls.diradjvls{ntind}=sort_times(avls.dirvls{ntind},avls.dirtimes');
%         end
        end
        
     