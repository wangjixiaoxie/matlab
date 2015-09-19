
function [ax]=plotdynamicscatter(sumdyn,norm,addx)

effratecomb=[];
pctratecomb=[];


col={ 'k' 'r'}
for ii=1:length(sumdyn)
%    colvl=mod(ii,4)
   if(ii==3)
       colvl=2
   else
       colvl=1;
   end
    smcr=sumdyn(ii);
    smcr.tms=smcr.tms-smcr.tms(1)+1;
    
      if(norm)
          plotpct=smcr.pct/smcr.pct(1);
      else
          plotpct=smcr.pct;
      end
    ind=find(smcr.contreff>0)
    
   hold on;
    plot(smcr.tms,plotpct,'Color',col{colvl})
    
    ind=find(smcr.contreff>0)
    clear tmdiff effrate pctrate
    
    if(length(ind)>=2)
        tms=smcr.tms(ind);
        pct=smcr.pct(ind);
        mneff=mean([smcr.targeff(ind);smcr.contreff(ind)]);
        mneffnorm=mneff/mneff(1);
        for kk=2:length(tms)
                tmdiff=tms(kk)-tms(kk-1);
                effrate(kk-1)=(mneffnorm(kk)-mneffnorm(kk-1))./tmdiff;
                pctrate(kk-1)=(pct(kk)-pct(kk-1))./tmdiff;
               
        end
    else
        effrate=[];
        pctrate=[];
    end
    effratecomb=[effratecomb effrate]
    pctratecomb=[pctratecomb pctrate]
    
end
figure
plot(effratecomb,pctratecomb,'k.')

     