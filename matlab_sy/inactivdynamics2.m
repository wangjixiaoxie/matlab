%made for plotting inactiv_dynamics
%inactiv_dynamics.m
function [ax]=inactivdynamics2(shsall, sumshs, datsource)

figure
subplot(211)
for ii=1:length(shsall)
    sl=shsall(ii)
    if (datsource=='sub')
        dat=sl.subruns
    else
        dat=sl.shiftruns
    end
   if(~isempty(dat)&length(dat)>1)
        ofst_tm=sl.flrtmvec-sl.flrtmvec(dat(1));
   end
        x=ofst_tm(dat);
   plot(x+1,sl.pct(dat),'k-');
   hold on;
end

plot([1:6], sumshs.pct,'c','Linewidth',2);
plot([1:6;1:6],[sumshs.pct-sumshs.stderrpct;sumshs.pct+sumshs.stderrpct],'c');


subplot(212)
for ii=1:length(shsall)
    %calc mean eff
    
    sl=shsall(ii)
    
    nts=sl.allnote
    mneff=mean(sl.effvls(nts,:),1)
        
   if(~isempty(dat)&length(dat)>1)
        ofst_tm=sl.flrtmvec-sl.flrtmvec(dat(1));
   end
        x=ofst_tm(dat);
   plot(x+1,mneff(dat),'k-');
   hold on;
  

end
 plot([1:6], sumshs.meaneff,'c','Linewidth',2);
plot([1:6;1:6],[sumshs.meaneff-sumshs.stderreff;sumshs.meaneff+sumshs.stderreff],'c');