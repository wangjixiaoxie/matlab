%rewritten extensively 4.16.09,
%to incorporate new figure requirements.

%takes sumdyn,
%and a ps (plotstruct).
%ps.minx - minimum x value to include
%ps.maxx - maximum x value to include.
%ps.col,  col={ 'k' 'r','c'}
%ps.colvec
%ps.type, (pct or off or dis)
%ps.addx
%ps.excludebs
%ps.plotavg=1
%ps.comb
%ps.normeff=0
%ps.plot_type, 1 is for shift, 2 is for asymp, 3 is for rev.


function [ctinds,meanmu,meanac]=plotcombdynamics2_stim(sumdynin,sumbs,ps)
% tm_matrix=1:.2:ps.maxx
if(ps.plotsum)
    
    axes(ps.ax(1));
end
axis(ps.plotbnds);

  
    sdynln=length(sumdynin)
    bsindcomb=[];
    sumdyncomb=[];
    outstruct=[];
    ctinds=[];
   fixct=0;
    ct=1;
    ln=ps.maxx
    if(ps.roundtimes)
            tmdiff=1;
             
            outvlaczcomb=zeros(sdynln,ln+1);
            outvlmuzcomb=zeros(sdynln,ln+1);
        
    else
        
            tmdiff=.1;
            outvlaczcomb=zeros(sdynln,(ln+1)*tmdiff);
            outvlmuzcomb=zeros(sdynln,(ln+1)*tmdiff);
    end
    
    
    for ii=1:length(sumdynin) 
        smcr=sumdynin(ii);
        shiftind=find(smcr.exadjtms>0);
%         aczvl(basind)=0;
%         muzvl(basind)=0;
        crsumbs=sumbs(smcr.bsnum);
        
        if(~isempty(smcr.acz))
                
                    acpctvl=smcr.pct(shiftind);
                    mupctvl=100-smcr.pct(shiftind);
                    aczvl=smcr.acz(shiftind)
                    
         end
            if(ps.flip)
                if(smcr.drxn=='do')
                    if(ps.plot_type=='pct')
%                         acpctvl=-acpctvl
%                         mupctvl=-mupctl
                        aczvl=-aczvl;
                    end
                end
            end
       
        
       if(ps.roundtimes)
            %this averages values on the same day
            adjtms=ceil(smcr.exadjtms(shiftind));
            [adjtms,acpctvl,mupctvl,aczvl]=adj_vls(adjtms,acpctvl,mupctvl,aczvl)
            
            lncomb=ln;
            combvls=[0:tmdiff:lncomb]
        else
            adjtms=smcr.exadjtms(shiftind);
           
            lncomb=(ln)*tmdiff;
            combvls=[0:tmdiff:lncomb]
       end
        
       if(ps.plotfixedpoints)
           fixct=fixct+1;
           for crtmwin=1:length(ps.tmwins)
               crtms=ps.tmwins{crtmwin};
               [vls,matchinds]=intersect(adjtms,crtms);
               if(~isempty(matchinds))
                    meanmupct(fixct,crtmwin)=mean(calcmeanstder2(mupctvl(matchinds)));
                    meanacpct(fixct,crtmwin)=mean(calcmeanstder2(acpctvl(matchinds)));
                    meanacz(fixct,crtmwin)=mean(calcmeanstder2(aczvl(matchinds)))
              else
                  meanmupct(fixct,crtmwin)=NaN;
                  meanacpct(fixct,crtmwin)=NaN;
                  meanacz(fixct,crtmwin)=NaN;
              end
           end
            plotinds=ps.xvls;
           notnaind=find(~isnan(meanmupct(fixct,:)));
           
           %make sure that there are values at first and last timepint
           nomatchvls=find(~ismember([1 length(ps.tmwins)],notnaind));
           
           if(isempty(nomatchvls))
           ctinds=[ctinds fixct]
               axes(ps.axmot)
           plot([plotinds(notnaind)],[meanmupct(fixct,notnaind)],'Color','k','Marker','o','MarkerSize', 4);
            hold on;
            if(ps.plotlman)
                axes(ps.axlman);
                plot([plotinds(notnaind)],[meanacpct(fixct,notnaind)],'color','r','Marker','o','MarkerSize', 4);
                hold on;
            end
            if(ps.plotacz)
                axes(ps.axacz)
                plot([plotinds(notnaind)],[meanacz(fixct,notnaind)],'Color','k','Marker','o','MarkerSize', 4);
                hold on;
               
            end
           end
       end
      
       if(ps.interptozero)
            tm_matrix=[0:tmdiff:ps.maxx];
            startvl=1;
       else
           tm_matrix=[min(adjtms):tmdiff:ps.maxx]
           
       end
 if(ps.plotsum)      
if(~isempty(adjtms))
 if((max(adjtms)>=ps.minx)&length(adjtms)>1)
    [acvls]=interp1(adjtms,aczvl,tm_matrix);
    [muvls]=interp1(adjtms,muzvl,tm_matrix);
       
        startvl=tm_matrix(1)/tmdiff+1;
        outvlaczcomb(ct,startvl:length(combvls))=acvls
        outvlmuzcomb(ct,startvl:length(combvls))=muvls
        ct=ct+1;
            if(ps.plotraw)
                adjtms=smcr.exadjtms(shiftind);
%             plot([0:ln-1],outvlaczcomb(ct-1,:),'Color',ps.ac_col,'Linewidth',2)
                tmsind=find(adjtms<=ps.maxx);
                [out,sortind]=sort(adjtms(tmsind));
                tmsind=tmsind(sortind);
                if(ps.plotlman)
                plot(tm_matrix,acvls,'Color','k','Linewidth',2)
                end
                hold on;
         
%                 plot(tm_matrix,acvls,'o','MarkerSize',5,'MarkerFaceColor',ps.ac_col)
%                 plot(tm_matrix,muvls,'o','MarkerSize',5,'MarkerFaceColor',ps.mu_col)
%             hold on;
                plot(tm_matrix,muvls,'Color','k','Linewidth',2)
%             plot([0:ln-1],outvlmuzcomb(ct-1,:),'Color',ps.mu_col,'Linewidth',2)
             end
%     end
        
      end
end
 end
end
 
 [outmumn]=calcmeanstder2(meanmupct(ctinds,:));
 [outacmn]=calcmeanstder2(meanacpct(ctinds,:));
 [outaczmn]=calcmeanstder2(meanacz(ctinds,:));
 plotinds=ps.xvls;
 for kk=1:length(plotinds)
     axes(ps.axmot)
 plot([plotinds(kk)-.2 ;plotinds(kk)+.2],[outmumn(kk);outmumn(kk)],'Color',ps.mu_col,'Linewidth',4);
  axes(ps.axlman)
 plotinds=ps.xvls;
 if(ps.plotlman)
 plot([plotinds(kk)-.2;plotinds(kk)+.2],[outacmn(kk);outacmn(kk)],'Color','r','Linewidth', 4);
 end
 end
 
 

 inds=[1:length(combvls)]


function [normout]=normvls(vls,norm)
    if(norm)
        normout=vls./vls(1);
    else
        normout=vls;
    end

function [yvl,yvl2,equalflag]=geteffvl(smcr,indx,ps);
     if(ps.comb&ps.normeff)
          %first normalize each
          normtargeff=norm_onevls(smcr.targeff(indx));
          normctrleff=norm_onevls(smcr.contreff(indx));
          yvl=mean([normtargeff;normctrleff]);
          yvl2=yvl;
          equalflag=1;
     elseif(~ps.comb&~ps.normeff)
            yvl=smcr.targeff(indx);
            yvl2=smcr.contreff(indx);
            equalflag=0;
     elseif(ps.comb&~ps.normeff)
            yvl=mean([smcr.targeff(indx); smcr.contreff(indx)]);
            yvl2=yvl;
            equalflag=1;
     elseif(~ps.comb&ps.normeff)
             yvl=norm_onevls(smcr.targeff(indx));
             yvl2=norm_onevls(smcr.contreff(indx));
             equalflag=0;
     end
                
function [yout]=norm_onevls(yin);
    yout=yin/mean(yin);
    
    
    function [slope]= calcslope(xin,yin);
        
        s=polyfit(xin,yin',1);
        slope=s(1);
     
 %find matching adjtms
 %say that adjtms input is 12334556
 function[out_tms,outacpct,outmupct,outacz]=adj_vls(adjtms,acpct,mupct,inz)
     
     %uniquetms=123456
     %ind=123568
     [uniquetms,ind]=unique(adjtms);
     for ii=1:length(uniquetms)
         ind=find(adjtms==uniquetms(ii))
         out_tms(ii)=uniquetms(ii);
         outacpct(ii)=mean(acpct(ind));
         outmupct(ii)=mean(mupct(ind));
         outacz(ii)=mean(inz(ind));
         
      end
     
    
       
      
     
     