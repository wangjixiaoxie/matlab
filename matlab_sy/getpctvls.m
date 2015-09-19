
%this is done for one set of shiftruns at a time.
function [ypct,ypcter]=getpctvls(avls,ntind,shs); 
    %do pct for all the runs
    subrns=shs.subruns
    acz=avls.acz(ntind,:)
    muz=avls.muz(ntind,:)
    yac=acz(subrns);
    ymu=muz(subrns);
    ydiff=-shs.acz(subrns)+shs.muz(subrns);
        ydiffer=(shs.acerrz(subrns)+shs.muerrz(subrns))/2;
        if(shs.mxshift)
            ydiff=-ydiff;
           
        end
        ypct=ydiff./yac*100;
       
        ypcter=ypct.*(ydiffer./(ymu +ydiffer./yac))/2
        
        
       