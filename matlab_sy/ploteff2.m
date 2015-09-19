%11.16.08
%this is called by shiftsumplot2, shiftplot2
%axin has length of number of subruns
function [axvec]=ploteff2(shs,axin)
      
      
      
      
        for ii=1:length(shs.subruns)
            if(~isempty(shs.subruns{ii}))
                axes(axin(ii))
                subrns=shs.subruns{ii}
                flrtmvc=shs.flrtmvec;
                ofst_tm=flrtmvc-flrtmvc(subrns(1));
                x=ofst_tm(subrns);
                
                    yeff=shs.effvls([shs.allnote],subrns);
                
            if(shs.ntype=='targ')    
                plot(x,yeff,'k-')
                hold on;
            end   
            box off;
        
            end
        end
       
     