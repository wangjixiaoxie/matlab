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

%NEED TO CREATE AN OUTPUT

function [outstruct]=plotdynamics3(sumdyn,ps)
cnt=0;
bsindcomb=[];
sumdyncomb=[];
if (~isfield(ps,'excludebs'))
    ps.excludeind=11;
else
    ps.excludeind=ps.excludebs;
end
    

% first find the number of traces
for ii=1:length(sumdyn)
    if((max(sumdyn(ii).tms))>=ps.minx&ii~=ps.excludeind)
        bsindcomb=[bsindcomb sumdyn(ii).bsnum];
        sumdyncomb=[sumdyncomb ii];
        cnt=cnt+1;
    end
end


outvlcomb=zeros(cnt,ps.maxx);
outvleffcomb=zeros(cnt,ps.maxx);
cnt=1;
for ii=1:length(sumdyn)
    crcol=ps.col{ps.colvec(ii)}
    smcr=sumdyn(ii);
    if(ps.type=='off')
        inyvl=normvls(smcr.normoff,ps.norm);
    elseif(ps.type=='pct')
        inyvl=normvls(smcr.pct,ps.norm);
    elseif(ps.type=='dis')
        inyvl=normvls(smcr.asympdist,ps.norm);
    end
    smcr.tms=smcr.tms-smcr.tms(1)+1;
    if(exist('ps.addx'))
        smcr.tms=smcr.tms+addx;
    end

    if((max(smcr.tms)>=ps.minx)&(ii~=ps.excludeind))
        indx=find(smcr.tms<=ps.maxx);
        for plotnum=1:2
            subplot(1,2,plotnum)
            if(plotnum==1)
                yvl=inyvl;
                equalflag=1;
            else
                [yvl,yvl2,equalflag]=geteffvl(smcr,indx,ps);
                yvl=1./yvl;
                yvl2=1./yvl2;
            end
            if(equalflag==0)
                numyvl=2
            else
                numyvl=1;
            end
            for jj=1:numyvl
                if(jj==2)
                    yvl=yvl2;
                end
                [outvltms,outvls]=interpvls3(smcr.tms(indx),yvl(indx));
                plot(smcr.tms(indx),outvls(smcr.tms(indx)),'.','Color',crcol)
                hold on;
                plot(outvltms,outvls,'Color',crcol)
                slope=calcslope(smcr.tms(indx),yvl(indx));
                equalflag=1;
            end
            if(ps.plotavg) 
                ln=length(outvls)
                diff=ps.maxx-ln;
                if(diff)
                    outvls(ln+1:ln+diff)=NaN
                end
                if plotnum==1
                    outvlcomb(cnt,:)=outvls(1:ps.maxx);
                    yslope(cnt)=slope
                else
                    outvleffcomb(cnt,:)=outvls(1:ps.maxx);
                    effslope(cnt)=slope;
                    cnt=cnt+1;
                end
            end
            
            
            
        end
    end
end
            
        
     
%   subplot(122)    
%     if(~isempty(indx))
%         %rewriting this so that there are two options
%         %1. plot control and target separately or joined
%         
%         %2. if separate, normalization is around 1...if together,
%         %normalization is around 1.  (no normalization means straight
%         %mean??)
%         for ii=[ps.plot]
%             if
%         
%         if(ps.eff=='both')
%             mneff=mean([smcr.targeff(indx);smcr.contreff(indx)]);
%         elseif(ps.eff=='cont')
%             mneff=smcr.contreff(indx);
%         else
%             mneff=smcr.targeff(indx);
%         end
%         normeff=normvls(mneff,ps.norm);
%         [outvltms,outvleff]=interpvls3(smcr.tms(indx),normeff(indx));
%         plot(smcr.tms(indx),outvleff(smcr.tms(indx)),'.','Color',crcol);
%         hold on;
%         plot(outvltms,outvleff,'Color',crcol);
%         if(ps.plotavg)
%            ln=length(outvleff)
%            diff=ps.maxx-ln;
%            if(diff)
%                outvleff(ln+1:ln+diff)=NaN
%            end
%           
%             
%             outvleffcomb(cnt,:)=outvleff(1:ps.maxx);
%             cnt=cnt+1;
%         end
%      end
%    end
%    
% end

%incredibly inefficient ways of calculating mean values

[outeffmn,outeffstd]=calcmeanstder(outvleffcomb);
[outpctmn,outpctstd]=calcmeanstder(outvlcomb);

subplot(121)
plot([1:ps.maxx],outpctmn,'Color','g','Linewidth',4);
plot([1:ps.maxx;1:ps.maxx],[outpctmn-outpctstd;outpctmn+outpctstd],'g')
subplot(122)
plot([1:ps.maxx],outeffmn,'Color','g','Linewidth',4);
plot([1:ps.maxx;1:ps.maxx],[outeffmn-outeffstd;outeffmn+outeffstd],'g')


outstruct.effvls=outvleffcomb;
outstruct.yvls=outvlcomb;
outstruct.effmn=outeffmn;
outstruct.effsd=outeffstd;
outstruct.yvlmn=outpctmn;
outstruct.yvlsd=outpctstd;
outstruct.bsvls=bsindcomb;
outstruct.dynind=sumdyncomb;
outstruct.yslope=yslope;
outstruct.effslope=effslope;

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