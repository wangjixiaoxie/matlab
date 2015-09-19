%plotalldays_inactivation.m
%written 8.20 in order to show each day's raw data which goes into overall
%stats.
% function [ax]=inactivplot3c(bs,bsind,rawvl,ppf,contrl)
function [ax]=shiftsumplot(shft_s,ppf,targ)
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

figure
for ii=1:length(shft_s(indtoplot))
    shs=shft_s(indtoplot(ii));
    subind=shs.subind;
      subx=shs.offset(subind)-shs.offset(subind(1))+1;
     plot(subx,shs.pctrcv(subind),'k');
     box off;
     hold on;
     plot([subx';subx'],[shs.pctrcv(subind)+shs.pcterrz(subind);shs.pctrcv(subind)-shs.pcterrz(subind)],'k');
 end
     
   