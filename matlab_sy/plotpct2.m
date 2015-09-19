function [axvec]=plotpct2(shs, axvec)
  
      for ii=1:length(shs.subruns)
        if ~isempty(shs.subruns{ii})
            axes(axvec(ii));
        
            subrns=shs.subruns{ii}
            flrtmvc=shs.flrtmvec;
        
            ofst_tm=flrtmvc-flrtmvc(subrns(1));
            x=ofst_tm(subrns);
            yac=shs.acz(subrns);
            ymuc=shs.muz(subrns);
            yacer=shs.acerrz(subrns);
            ydiff=-shs.acz(subrns)+shs.muz(subrns);
            ydiffer=(shs.acerrz(subrns)+shs.muerrz(subrns))/2;
            if(shs.mxshift(ii))
                ydiff=-ydiff;
           
            end
            ypct=ydiff./yac*100;
       
            ypcter=ypct.*(ydiffer./(ymuc +ydiffer./yac))/2
            if(shs.ntype=='targ')    
                plot(x,ypct,'k-')
                hold on;
                plot([x'; x'],[(ypct+ypcter); (ypct-ypcter)],'k');
            end   
            box off;
        end
        
      end
