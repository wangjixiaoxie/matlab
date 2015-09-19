
function [ax]=plotdynamics(sumdyn,norm,numx,plotavg,offz,colind,addx)


    outvlpctcomb=[];
    outvleffcomb=[];
    

col={ 'k' 'r','c'}
for ii=1:length(sumdyn)
%    colvl=mod(ii,4)
   if(ii==3)
       colvl=2
   elseif (ii==13||ii==14)
       colvl=3;
   else
       colvl=1;
   end
    smcr=sumdyn(ii);
    if(offz)
        if(norm)
            yvl=smcr.normoff/smcr.normoff(1);
        else
            yvl=smcr.normoff;
        end
    else
        if(norm)
          yvl=smcr.pct/smcr.pct(1);
        else
          yvl=smcr.pct;
        end
    end
    smcr.tms=smcr.tms-smcr.tms(1)+1;
    if(exist('addx'))
        smcr.tms=smcr.tms+addx;
    end
    subplot(211)
      
     
      
   
   hold on;
   if((length(smcr.tms)>1)&(max(smcr.tms)>2)) 
      
       [outvltms,outvlpct]=interpvls2(smcr.tms,yvl);
       indx=find(smcr.tms<=numx); 
          plot(smcr.tms(indx),yvl(indx),'.','Color',col{colvl})
          hold on;
        if(length(outvltms)>=numx)
            plot(outvltms(1:numx),outvlpct(1:numx),'Color',col{colvl})
        end
            if(plotavg)
            [outvltms,outvlpct]=interpvls2(smcr.tms,yvl);
            if(length(outvltms)>=numx)
                outvlpctcomb=[outvlpctcomb;outvlpct(1:numx)];
            end
        end
     
  subplot(212)
%    plot(smcr.adjshiftms,mean(smcr.logtargeff,'k');
   hold on
  
    
    ind=find(smcr.contreff>0)
    if(~isempty(ind))
      mneff=mean([smcr.targeff(ind);smcr.contreff(ind)]);
      if(norm)
          plotmneff=mneff/mneff(1);
      else
          plotmneff=mneff;
      end
         [outvltms,outvleff]=interpvls2(smcr.tms,plotmneff);

    plot(smcr.tms(indx),plotmneff(indx),'.','Color',col{colvl});
    hold on;
    if(length(outvltms)>=numx)
        plot(outvltms(1:numx),outvleff(1:numx),'Color',col{colvl});
    end
   else
    mneff=smcr.targeff;
    if(norm)
        plotmneff=mneff/mneff(1);
    else
        plotmneff=mneff;
%     end
%     plot(smcr.tms(indx),plotmneff ,'r');
    end
    end
  
    if(plotavg)
%         [outvltms,outvleff]=interpvls2(smcr.tms,plotmneff);
        if(length(outvltms)>=numx)
            outvleffcomb=[outvleffcomb;outvleff(1:numx)];
                end
     end
   end
   
end

outeffmn=mean(outvleffcomb);
outeffstd=std(outvleffcomb)/sqrt(length(outvleffcomb(:,1)));
outpctmn=mean(outvlpctcomb);
outpctstd=std(outvlpctcomb)/sqrt(length(outvlpctcomb(:,1)));

subplot(211)
plot([1:numx],outpctmn,'k','Linewidth',2);
plot([1:numx;1:numx],[outpctmn-outpctstd;outpctmn+outpctstd],'k')
subplot(212)
plot([1:numx],outeffmn,'k','Linewidth',2);
plot([1:numx;1:numx],[outeffmn-outeffstd;outeffmn+outeffstd],'k')


