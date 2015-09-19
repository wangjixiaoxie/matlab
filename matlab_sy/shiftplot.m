%plotalldays_inactivation.m
%written 8.20 in order to show each day's raw data which goes into overall
%stats.
% function [ax]=inactivplot3c(bs,bsind,rawvl,ppf,contrl)
function [ax]=shiftplot(shft_s,ppf,targ)
    %in first run, setting this up as a three column subplout
    %first column is raw zscore values for acsf/muscimol
    %second column only shows points beginning with first point which is
    %within one standard dev of maximum shift
    %2nd column is amount of actual pitchshift
    %3rd column is the percent recovery.
 figcount=1;
 ncol=4;
 
 %find which values are targeted and which values are controls.
 
 for ii=1:length(shft_s)
     if(shft_s(ii).ntype=='targ')
         targvec(ii)=1
     else
         targvec(ii)=0;
     end
 end
     
     
if (targ=='cont')
    indtoplot=find(targvec==0);
else
    indtoplot=find(targvec==1)

end


for ii=1:length(shft_s(indtoplot))
     shs=shft_s(indtoplot(ii));
     subind=shs.subind;
     pn=mod(ii,ppf);
     if ~pn
         pn=ppf;
     end
     if (pn==1)
         figure;
         cnt=1;
     end
     %each ii is offset by ncols for total count.
     axcnt=ncol*(ii-1)+1;
     pltind=mod(axcnt,ppf*ncol);
     %plot raw data
     ax(axcnt)=subplot(ppf,4,pltind);
     plot(shs.offsetall,shs.aczall,'k-')
     hold on;
     plot(shs.offsetall,shs.muzall,'r-');
     plot([shs.offsetall'; shs.offsetall'],[(shs.aczall+shs.acerrzall); (shs.aczall-shs.acerrzall)],'k');
     plot([shs.offsetall'; shs.offsetall'],[(shs.muzall+shs.muerrzall); (shs.muzall-shs.muerrzall)],'r');
     box off;
     %plotsubdata
     subx=shs.offset(subind)-shs.offset(subind(1))+1;
     subac=shs.acz(subind);
     submu=shs.muz(subind);
     subacerr=shs.acerrz(subind);
     submuerr=shs.muerrz(subind);
     box off;
     title(shs.bname);
     ax(axcnt+1)=subplot(ppf,4,pltind+1);
     plot(subx, subac,'k-')
     hold on;
     plot(subx, submu,'r-');
     plot([subx'; subx'],[(subac+subacerr); (subac-subacerr)],'k');
     plot([subx'; subx'],[(submu+submuerr); (submu-submuerr)],'r');
     box off;
     %plot offset
     ax(axcnt+2)=subplot(ppf,4,pltind+2);
     if(shs.dir=='up')
        plot(subx,shs.diffz(subind),'k')
        hold on;
        plot([subx';subx'],[shs.diffz(subind)+shs.differrz(subind);shs.diffz(subind)-shs.differrz(subind)],'k')
     else
        plot(subx,-shs.diffz(subind),'k')
        hold on;
         plot([subx';subx'],[-shs.diffz(subind)+shs.differrz(subind);-shs.diffz(subind)-shs.differrz(subind)],'k')
     end
    
     box off;
     %plot percent
     ax(axcnt+3)=subplot(ppf,4,pltind+3);
     plot(subx,shs.pctrcv(subind),'k');
     box off;
     hold on;
     plot([subx';subx'],[shs.pctrcv(subind)+shs.pcterrz(subind);shs.pctrcv(subind)-shs.pcterrz(subind)],'k');
 end
     
   