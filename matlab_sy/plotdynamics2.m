%rewritten extensively 4.16.09,
%to incorporate new figure requirements.

%takes sumdyn,
%and a ps (plotstruct).
%ps.norm - whether to normalize or not
%ps.minx - minimum x value to include
%ps.maxx - maximum x value to include.
%ps.colors,  col={ 'k' 'r','c'}
%ps.colvec
%ps.type, (pct or off)
%ps.addx

%NEED TO CREATE AN OUTPUT

function [outstruct]=plotdynamics2(sumdyn,ps)
cnt=0;
bsindcomb=[];
sumdyncomb=[];
if (~isfield(ps,'excludebs'))
    ps.excludeind=11;
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
    else
        inyvl=normvls(smcr.pct,ps.norm);
    end
    smcr.tms=smcr.tms-smcr.tms(1)+1;
    if(exist('ps.addx'))
        smcr.tms=smcr.tms+addx;
    end
    subplot(121)
    hold on;
if((max(smcr.tms)>=ps.minx)&(ii~=ps.excludeind))
     indx=find(smcr.tms<=ps.maxx);
     if(~isempty(indx))
        [outvltms,outvls]=interpvls3(smcr.tms(indx),inyvl(indx));
        plot(smcr.tms(indx),outvls(smcr.tms(indx)),'.','Color',crcol)
        hold on;
        plot(outvltms,outvls,'Color',crcol)
        
       if(ps.plotavg) 
           ln=length(outvls)
           diff=ps.maxx-ln;
           if(diff)
               outvls(ln+1:ln+diff)=NaN
           end
           outvlcomb(cnt,:)=outvls(1:ps.maxx);
       end
     end
     
  subplot(122)    
    if(~isempty(indx))
        if(ps.eff=='both')
            mneff=mean([smcr.targeff(indx);smcr.contreff(indx)]);
        elseif(ps.eff=='cont')
            mneff=smcr.contreff(indx);
        else
            mneff=smcr.targeff(indx);
        end
        normeff=normvls(mneff,ps.norm);
        [outvltms,outvleff]=interpvls3(smcr.tms(indx),normeff(indx));
        plot(smcr.tms(indx),outvleff(smcr.tms(indx)),'.','Color',crcol);
        hold on;
        plot(outvltms,outvleff,'Color',crcol);
        if(ps.plotavg)
           ln=length(outvleff)
           diff=ps.maxx-ln;
           if(diff)
               outvleff(ln+1:ln+diff)=NaN
           end
          
            
            outvleffcomb(cnt,:)=outvleff(1:ps.maxx);
            cnt=cnt+1;
        end
     end
   end
   
end

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

function [normout]=normvls(vls,norm)
    if(norm)
        normout=vls./vls(1);
    else
        normout=vls;
    end
