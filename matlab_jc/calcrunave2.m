 function [outmean,outstd,out_tms]=calcrunave2(vlsout,WS)
    
     WSint=floor(WS/2);
%      PTint=floor(PWS/2);
     stind=WSint+1;
%      ptint=PTint+1;
    for ii=1:length(vlsout)
        crvls=vlsout{ii}
        ln=length(crvls(:,1))
        vec=WSint+1:ln-WSint;
        for jj=1:length(vec)
             vlind=vec(jj);
             vecind=vlind-WSint:vlind+WSint
             out_tms{ii}(jj)=mean(crvls(vecind,1));
             outmean{ii}(jj)=mean(crvls(vecind,2));
             outstd{ii}(jj)=std(crvls(vecind,2));
             %could add something to show what mean threshold by separating
             %hits and misses
        end
   
    end
        
